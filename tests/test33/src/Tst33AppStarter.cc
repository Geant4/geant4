//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: Tst33AppStarter.cc,v 1.12 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33AppStarter.cc
//
// ----------------------------------------------------------------------

#include "Tst33RunAction.hh"
#include "Tst33AppStarter.hh"
#include "Tst33PhysicsList.hh"
#include "Tst33PrimaryGeneratorAction.hh"
#include "Tst33DetectorConstruction.hh"
#include "G4UImanager.hh"
#include "G4GeometrySampler.hh"
#include "Tst33VGeometry.hh"
#include "Tst33SlobedConcreteShield.hh"
#include "Tst33ConcreteShield.hh"
#include "Tst33ParallelGeometry.hh"
#include "G4CellScorerStore.hh"
//#include "G4CellStoreScorer.hh"
#include "Tst33IStoreBuilder.hh"
#include "Tst33WeightWindowStoreBuilder.hh"
#include "Tst33ScorerBuilder.hh"
#include "Tst33VEventAction.hh"
#include "Tst33VisRunAction.hh"
#include "Tst33VisEventAction.hh"
#include "Tst33TimedEventAction.hh"
#include "G4VIStore.hh"
#include "G4CellScorer.hh"
//#include "G4ScoreTable.hh"
#include "Tst33ScoreTable.hh"
#include "Tst33TimedApplication.hh"
#include "Tst33VisApplication.hh"
#include "G4ProcessPlacer.hh"
#include "Tst33WeightChangeProcess.hh"
#include "G4PlaceOfAction.hh"
#include "G4WeightWindowAlgorithm.hh"
#include "G4VWeightWindowStore.hh"

// customized scorer and store
#include "Tst33CellScorer.hh"
#include "Tst33CellScorerStore.hh"

// a score table
//#include "G4ScoreTable.hh"
#include "Tst33ScoreTable.hh"



Tst33AppStarter::Tst33AppStarter()
  : 
  fMessenger(*this),
  fMassGeometry(0),
  fSampleGeometry(0),
  fParallelGeometry(0),
  fApp(0),
  fDetectorConstruction((new Tst33DetectorConstruction)),
  fSampler(0),
  fScorer(0),
  fCell_19_Scorer(0),
  fIStore(0),
  fWWStore(0),
  fConfigured(false),
  fWeightroulette(false),
  fTime(0),
  fChangeWeightPlacer(0),
  fWeightChangeProcess(0),
  fWWAlg(0),
  parallel_geometry(false),
  forceCoupled(false)
{
  if (!fDetectorConstruction) {
    NewFailed("Tst33AppStarter", "Tst33DetectorConstruction");
  }
  physlist = new Tst33PhysicsList;
  fRunManager.SetUserInitialization(physlist); 
  fRunManager.SetUserAction(new Tst33PrimaryGeneratorAction);  
  fRunManager.SetUserInitialization(fDetectorConstruction);
  fRunManager.SetUserAction(new Tst33RunAction);
}
Tst33AppStarter::~Tst33AppStarter()
{
  if (fSampler) {
    delete fSampler;
  }
  if (fParallelGeometry) {
    delete fParallelGeometry;
  }
  if (fWeightroulette) {
    fChangeWeightPlacer->RemoveProcess(fWeightChangeProcess);
    delete fWeightChangeProcess;
    fWeightChangeProcess = 0;
    delete fChangeWeightPlacer;
    fChangeWeightPlacer = 0;
  }
  if (fWWAlg) {
    delete fWWAlg;
  }
}

void Tst33AppStarter::CreateMassGeometry() {
  G4cout << " Creating Mass Geometry " << G4endl;
  if (fMassGeometry) {
    G4cout << "Tst33AppStarter::CreateMassGeometry: a geometry already exists!" << G4endl;
  } 
  else {
    fMassGeometry = new Tst33SlobedConcreteShield;
    if (!fMassGeometry) {
      NewFailed("CreateMassGeometry", "Tst33SlobedConcreteShield");
    }

    fMassGeometry->Construct();
    
    fDetectorConstruction->SetWorldVolume(&fMassGeometry->
					  GetWorldVolumeAddress());
  }

  parallel_geometry = false;

  if(forceCoupled) physlist->UseCoupledTransportation(true);

}

void Tst33AppStarter::CreateParallelGeometry() {
  if (fMassGeometry) {
    G4cout << "Tst33AppStarter::CreateParallelGeometry: a eometry already exists!" << G4endl;
  } 
  else {
    fMassGeometry = new Tst33ConcreteShield();
    if (!fMassGeometry) {
      NewFailed("CreateParallelGeometry", "Tst33ConcreteShield");
    }    
    
    fMassGeometry->Construct();
    
    fDetectorConstruction->SetWorldVolume(&fMassGeometry->
					  GetWorldVolumeAddress());
    G4String parallelName("ParallelBiasingWorld"); //ASO
    fParallelGeometry = new Tst33ParallelGeometry(parallelName,
						  &fMassGeometry->GetWorldVolumeAddress());
    if (!fParallelGeometry) {
      NewFailed("CreateParallelGeometry", "Tst33ParallelGeometry");
    }

    //    fParallelGeometry->Construct();

    fDetectorConstruction->RegisterParallelWorld(fParallelGeometry);
    //    fDetectorConstruction->RegisterParallelWorld(fParallelGeometry->GetParallelWorld());
//     fSampleGeometry = fParallelGeometry;
// //     G4GeometrySampler fSampler(&fParallelGeometry->GetWorldVolume(),
// // 			       "neutron");
// //    fSampler = new G4GeometrySampler(&fParallelGeometry->GetWorldVolume(),
// //			       "neutron");
//     fSampler = new G4GeometrySampler(&fParallelGeometry->GetWorldVolumeAddress(),
// 				     "neutron");

//     fSampler->SetParallel(true);

    physlist->AddParallelWorldName(parallelName);

  }

  parallel_geometry = true;

}

G4bool Tst33AppStarter::IsGeo_n_App(){
  G4bool geo_n_app = true;
  if (!fMassGeometry) {
    G4cout << "Tst33AppStarter::IsGeo_n_App: no geometry!" << G4endl;
    geo_n_app = false;
  }
  if (!fApp) {
    G4cout << "Tst33AppStarter::IsGeo_n_App: no application!" << G4endl;
    geo_n_app = false;
  }
  return geo_n_app;
}

G4bool Tst33AppStarter::CheckCreateScorer() {
  G4bool createscorer = true;
  if (!IsGeo_n_App()) {
    G4cout << "Tst33AppStarter::CheckCreateScorer: !IsGeo_n_App()!" << G4endl;
    createscorer = false;
  }
  if (fScorer) {
    G4cout << "Tst33AppStarter::CheckCreateScorer: scorer already exists!" << G4endl;
    createscorer = false;
  }
  if (fScorerStore) {
    G4cout << "Tst33AppStarter::CheckCreateScorer: scorer already exists!" << G4endl;
    createscorer = false;
  }
  return createscorer;
}

void Tst33AppStarter::CreateScorer() {
  // create a customized cell scorer store
  Tst33CellScorerStore tst33store;

  if (CheckCreateScorer()) { 
    Tst33ScorerBuilder sb;
    fScorerStore = sb.CreateScorer(fSampleGeometry,
 				   &fCell_19_Scorer);
//     fScorer = new G4CellStoreScorer(*fScorerStore);
//     if (!fScorer) {
//       NewFailed("CreateScorer","G4CellStoreScorer");
//     }
//     fSampler->PrepareScoring(fScorer);
//x    fEventAction->SpecialCellScorer(fCell_19_Scorer);
  }
}

void Tst33AppStarter::DeleteScorers() {
  //  delete fScorer;
  //  fScorer = 0;
  fScorerStore->DeleteAllScorers();
  fCell_19_Scorer = 0;
}

G4bool Tst33AppStarter::CheckCreateIStore() {
  G4bool createistore = true;
  if (!IsGeo_n_App()) {
    G4cout << "Tst33AppStarter::CheckCreateIStore: !IsGeo_n_App()!" << G4endl;
    createistore = false;
  }  
  if (fIStore) {
    G4cout << "Tst33AppStarter::CheckCreateIStore: an importance store already exists!" << G4endl;
    createistore = false;
  }
  return createistore;
}

G4bool Tst33AppStarter::CheckCreateWeightWindowStore() {
  G4bool createWWstore = true;
  if (!IsGeo_n_App()) {
    G4cout << "Tst33AppStarter::CheckCreateWeightWindowStore: !IsGeo_n_App()!" << G4endl;
    createWWstore = false;
  }  
  if (fWWStore) {
    G4cout << "Tst33AppStarter::CheckCreateWeightWindowStore: an weight window store already exists!" << G4endl;
    createWWstore = false;
  }
  return createWWstore;
}

void Tst33AppStarter::CreateIStore() {
  if (CheckCreateIStore()) {
    Tst33IStoreBuilder ib;
    fIStore = ib.CreateIStore(fSampleGeometry);
    fSampler->PrepareImportanceSampling(fIStore);
  }
}

void Tst33AppStarter::CreateWeightWindowStore(G4PlaceOfAction poa,
					      G4bool zeroWindow) {
  G4cout << " Weight window store creator 1 " << G4endl;
  if (CheckCreateWeightWindowStore()) {
    
    Tst33WeightWindowStoreBuilder wb;
    G4cout << " Weight window store creator 2 " << G4endl;
    fWWStore = wb.CreateWeightWindowStore(fSampleGeometry);
    G4cout << " Weight window store creator 3 " << G4endl;
    if (zeroWindow) {
      fWWAlg = new G4WeightWindowAlgorithm(1,1,100);
    }
    else {
      fWWAlg = new G4WeightWindowAlgorithm(5,3,5);
    }
    G4cout << " Weight window store creator 4 " << poa << G4endl;
    if(!fSampler) G4cout << " fSampler doesn't exist - why? " << G4endl; 
    fSampler->PrepareWeightWindow(fWWStore,
				  fWWAlg,
				  poa);
    //				  onBoundary);
    G4cout << " Weight window store creator 5 " << G4endl;



  }
}

G4bool Tst33AppStarter::CheckCreateWeightRoulette() {
  G4bool createWR = true;
  if (!IsGeo_n_App()) {
    G4cout << "Tst33AppStarter::CheckCreateWeightRoulette: !IsGeo_n_App()!" << G4endl;
    createWR = false;
  }  
  if (fWeightroulette) {
    G4cout << "Tst33AppStarter::CheckCreateWeightRoulette: weigyht roulette already exists!" << G4endl;
    createWR = false;
  }
  return createWR;
}

void Tst33AppStarter::AddWeightChanger(){
  
  fWeightChangeProcess = new Tst33WeightChangeProcess;
  fChangeWeightPlacer = new G4ProcessPlacer("neutron");
  fChangeWeightPlacer->AddProcessAsLastDoIt(fWeightChangeProcess);
}

void Tst33AppStarter::CreateWeightRoulette(G4int mode) {
  if (CheckCreateWeightRoulette()) {
    fWeightroulette = true;
    AddWeightChanger();
    if (mode==1) {
      // make sure that no relative weights below 1 occur
      fSampler->PrepareWeightRoulett(1, 1, 1);
    }
    else {
      // apply default weight eoulette:  see G4VSampler.hh
      fSampler->PrepareWeightRoulett();
    }
  }
}


G4bool Tst33AppStarter::CheckCreateApp() {
  G4bool  createapp(true);
  if (!fMassGeometry) {
    G4cout << "Tst33AppStarter::CheckApp: geometry not defined!" << G4endl;
    createapp = false;
  }
  if (fApp) {
    G4cout << "Tst33AppStarter::CheckApp: an application alredy exists!" << G4endl;
    createapp = false;
  }   
  return createapp;
}

void Tst33AppStarter::CreateVisApplication(){
  if (CheckCreateApp()) {
    fApp = new Tst33VisApplication;
    if (!fApp) {
      NewFailed("CreateVisApplication", "Tst33VisApplication");
    }
    fRunAction = fApp->CreateRunAction();
    fRunManager.SetUserAction(fRunAction);
    fEventAction = fApp->CreateEventAction();
    fRunManager.SetUserAction(fEventAction);    
    fRunManager.Initialize();

    if(parallel_geometry) {
      fSampleGeometry = fParallelGeometry;
      //     G4GeometrySampler fSampler(&fParallelGeometry->GetWorldVolume(),
      // 			       "neutron");
      //    fSampler = new G4GeometrySampler(&fParallelGeometry->GetWorldVolume(),
      //			       "neutron");
      fSampler = new G4GeometrySampler(&fParallelGeometry->GetWorldVolumeAddress(),
				       "neutron");
      
      fSampler->SetParallel(true);
    } else {
      fSampleGeometry = fMassGeometry;
      //    G4GeometrySampler fSampler(&fMassGeometry->GetWorldVolume(),"neutron");
      fSampler = new G4GeometrySampler(&fMassGeometry->GetWorldVolumeAddress(),"neutron");    
    }    

    G4UImanager::GetUIpointer()->ApplyCommand("/control/execute vis.mac");
  }
}

void Tst33AppStarter::CreateTimedApplication(G4int time) {
  if (CheckCreateApp()) {
    fTime = time;
    fApp = new Tst33TimedApplication(fTime);
    if (!fApp) {
      NewFailed("CreateTimedApplication","Tst33TimedApplication");
    }
    fEventAction = fApp->CreateEventAction();
    fRunManager.SetUserAction(fEventAction);    
    fRunManager.Initialize();

    if(parallel_geometry) {
      fSampleGeometry = fParallelGeometry;
      //     G4GeometrySampler fSampler(&fParallelGeometry->GetWorldVolume(),
      // 			       "neutron");
      //    fSampler = new G4GeometrySampler(&fParallelGeometry->GetWorldVolume(),
      //			       "neutron");
      fSampler = new G4GeometrySampler(&fParallelGeometry->GetWorldVolumeAddress(),
				       "neutron");
      
      fSampler->SetParallel(true);
    } else {
      fSampleGeometry = fMassGeometry;
      //    G4GeometrySampler fSampler(&fMassGeometry->GetWorldVolume(),"neutron");
      fSampler = new G4GeometrySampler(&fMassGeometry->GetWorldVolumeAddress(),"neutron");    
    }    
  }
}


void Tst33AppStarter::ClearSampling() {
  if (!fSampler) {
    G4cout << "Tst33AppStarter::ClearSampling: no sampler!" << G4endl;
  }
  else {
    fSampler->ClearSampling();
    if (fScorer) {
       DeleteScorers();
    }
    if (fIStore) {
      delete fIStore;
      fIStore = 0;
    }
    fEventAction->Clear();
    if (fWeightroulette) {
      fChangeWeightPlacer->RemoveProcess(fWeightChangeProcess);
      delete fWeightChangeProcess;
      fWeightChangeProcess = 0;
      delete fChangeWeightPlacer;
      fChangeWeightPlacer = 0;
    }
    if (fWeightChangeProcess) {
      fChangeWeightPlacer->RemoveProcess(fWeightChangeProcess);
      delete fWeightChangeProcess;
      fWeightChangeProcess = 0;
      delete fChangeWeightPlacer;
      fChangeWeightPlacer = 0;
    }
    if (fWWStore) {
      delete fWWStore;
      fWWStore = 0;
    }
    if (fWWAlg) {
      delete fWWAlg;
      fWWAlg = 0;
    }
    fWeightroulette = false;
    fConfigured = false;
  }
}

void Tst33AppStarter::ConfigureSampling(){
  if (!fConfigured) {
    if (!IsGeo_n_App()) {
      G4cout << "ConfigureSampling no geometry or application given!" << G4endl;
    }
    else {
      fConfigured = true;
      fSampler->Configure();
    }
  }
  else {
    G4cout << "Tst33AppStarter::ConfigureSampling(): already configured!" << G4endl;
  }
}

void Tst33AppStarter::Run(G4int nevents) {
  if (!fConfigured) {
    G4cout << "Not configured yet!" << G4endl;
  } 
  else {
    G4int n = nevents;
    if (fTime >0 ) {
      n = 10000000;
    }
    fRunManager.BeamOn(n);
  }	
}

void Tst33AppStarter::PostRun() {
  if (fConfigured) {
    if (fScorerStore) {
      //      G4ScoreTable sp(fIStore);
      //      B02ScoreTable sp(&aIstore);
      //xtest      Tst33ScoreTable sp(fIStore);
      //xtest      sp.Print(fScorerStore->GetMapGeometryCellCellScorer(), &G4cout);
    }
  }
  else {
    G4cout << "Tst33AppStarter::PostRun: not configured!" << G4endl;
  }
}


void Tst33AppStarter::NewFailed(const G4String &function, 
				const G4String &cl){
    G4String where="Tst33AppStarter::"+ function +"()";
    G4ExceptionDescription what;
    what << "new failed to create: "<< cl << "!";				
    G4Exception(where, "TST33-01",FatalException, what);
}





