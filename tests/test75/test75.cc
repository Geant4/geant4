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

#include "globals.hh"
#include "G4PhysicalConstants.hh"

#include "G4ios.hh"

#include "G4Material.hh"
#include "G4IsotopeVector.hh"
#include "G4ElementVector.hh"
#include "G4NistManager.hh"

#include "G4ProcessManager.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4HadronCrossSections.hh"
#include "G4VCrossSectionDataSet.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"

#include "G4Gamma.hh"
#include "G4ForceCondition.hh"

#include "G4CascadeInterface.hh"
#include "G4PhotoNuclearProcess.hh" // hadronic process
// #include "G4GammaNuclearReaction.hh" // hadronic interaction (aka "model")


#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"
#include "G4GRSSolid.hh"

#include "G4TrackingManager.hh"
#include "G4Region.hh"

#include "G4StateManager.hh"
#include "G4Navigator.hh"
#include "G4Timer.hh"

// random engine/seed settings
#include "Randomize.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TString.h"

#include "G4SystemOfUnits.hh"

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

int main(int argc, char** argv) {

  G4cout <<"========================================================" <<G4endl;
  G4cout <<"======      Gamma-Nuclear Test Start       ========" <<G4endl;
  G4cout <<"========================================================" <<G4endl;

  G4String  namePart = "gamma";
  G4String  nameMatRaw  = "Cu";
  ostringstream osMat(ios_base::out|ios_base::app);// to enable appending in output operations
  G4String  nameGen  = "Bertini";
  G4double  energy   = 0.;
  G4double  m_p      = 300.*MeV;
  G4int     nevt     = 1000;
  G4double theStep   = 0.01*micrometer; // should it not be 0. ?
  G4Material* material = 0;
  G4int  verbose  = 0; 

//
//
  G4long myseed   = 135799753;
  G4int  jobid = -1;
  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(myseed);

  // Track
  G4ThreeVector aPosition  = G4ThreeVector(0.,0.,0.);
  G4double      aTime      = 0. ;
  G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
  G4double nx = 0.0, ny = 0.0, nz = 0.0;

  // Control on input

  if (argc < 2) {
    G4cout << "Input file is not specified! Exit" << G4endl;
    exit(1);
  }

  std::ifstream* fin = new std::ifstream();
  G4String fname = argv[1];
  fin->open(fname.c_str());
  if( !fin->is_open()) {
    G4cout << "Input file <" << fname << "> does not exist! Exit" << G4endl;
    exit(1);
  }

  /*  G4UImanager* UI = G4UImanager::GetUIpointer();
      G4String command("/tracking/verbose 2");
      UI->ApplyCommand( command );
  */  

  // we need to be in the preinit state to create particles

  if(!G4StateManager::GetStateManager()->SetNewState(G4State_PreInit))
    G4cout << "G4StateManager PROBLEM! " << G4endl;


  // physics needs to be initialized before the 1st use of particle table,
  // because it constructs particles - otherwise the table is just empty
  //
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();  
  
  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();
  
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();
  
  G4BosonConstructor pBosonConstructor;
  pBosonConstructor.ConstructParticle();
  

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  // Geometry
  //
  G4Box* sFrame = new G4Box ("Box", 100.0*cm, 100.0*cm, 100.0*cm);
  G4LogicalVolume* lFrame = new G4LogicalVolume(sFrame,material,"Box",0,0,0);
  G4PVPlacement* pFrame = new G4PVPlacement(0,G4ThreeVector(),"Box",
                                            lFrame,0,false,0);
  
  assert(pFrame);

  G4Region* rFrame = new G4Region("Region"); // needed by tracking manager
  lFrame->SetRegion(rFrame);
  rFrame->AddRootLogicalVolume(lFrame);

  G4TrackingManager* trackManager = new G4TrackingManager;


  // ---- Read input file
  G4cout << "Available commands are: " << G4endl;
  G4cout << "#events" << G4endl;
  G4cout << "#particle" << G4endl;
  G4cout << "#energy(MeV)" << G4endl;
  G4cout << "#momentum(MeV/c)" << G4endl;
  G4cout << "#step(mm)" << G4endl;
  G4cout << "#material" << G4endl;
  G4cout << "#generator" << G4endl;
  G4cout << "#verbose" << G4endl;
  G4cout << "#position(mm)" << G4endl;
  G4cout << "#direction" << G4endl;
  G4cout << "#time(ns)" << G4endl; // why would I need this ?
//
// for parallel processing
//
  G4cout << "#randomSeed" << G4endl;
  G4cout << "#jobID" << G4endl;
//
  G4cout << "#run" << G4endl;
  G4cout << "#exit" << G4endl;

  G4String line, line1;
  G4bool end = true;
  for(G4int run=0; run<100; run++) {
    do {
      (*fin) >> line;
      G4cout << "Next line " << line << G4endl;
      if(line == "#particle") {
        (*fin) >> namePart;
      } else if(line == "#momentum(MeV/c)") {
        (*fin) >> m_p;
        m_p *= MeV;
      } else if(line == "#events") {
        (*fin) >> line1;
        istringstream is(line1);
        is >> nevt;
        G4cout << "nevt : " << nevt << G4endl;
      } else if(line == "#step(mm)") {
        (*fin) >> theStep;
        theStep *= mm;
      } else if(line == "#material") {
        (*fin) >> nameMatRaw;
      } else if(line == "#generator") {
        (*fin) >> nameGen;
        G4cout << nameGen << G4endl;
      } else if(line == "#run") {
        break;
      } else if(line == "#verbose") {
        (*fin) >> verbose;
        if (verbose > 0) {
          G4cout <<"test51"
                 <<" : verbose set to "
                 << verbose
                 << G4endl;
        }
      } else if(line == "#position(mm)") {
        (*fin) >> nx >> ny >> nz;
        aPosition = G4ThreeVector(nx*mm, ny*mm, nz*mm);
      } else if(line == "#direction") {
        (*fin) >> nx >> ny >> nz;
        if(nx*nx + ny*ny + nz*nz > 0.0) {
          aDirection = G4ThreeVector(nx, ny, nz);
          aDirection = aDirection.unit();
	}
      } else if(line == "#time(ns)") {
        (*fin) >> aTime;
        aTime *= ns;
      } 
//
// needed for parallel processing
//
      else if ( line == "#randomSeed" )
      {
        (*fin) >> myseed;
	CLHEP::HepRandom::setTheSeed(myseed);
	G4cout << "###### Set Random Seed to " << myseed << "     #####" << G4endl;      
      }
      else if ( line == "#jobID" )
      {
        (*fin) >> jobid ;
      }
//
      else if(line == "#exit") 
      {
        end = false;
        break;
      }
    } while(end);

    if(!end) break;

    G4cout << "###### Start new run # " << run << "     #####" << G4endl;

    osMat.clear();
    osMat.str("G4_");
    osMat << nameMatRaw;
    G4String nameMat = osMat.str();

    G4cout << "###### Material: " << nameMat << " derived from " << nameMatRaw << G4endl;

    material = G4NistManager::Instance()->FindOrBuildMaterial(nameMat);

    if (!material) {
      G4cout << "Material <" << nameMat << "> is not found" << G4endl;
      exit(1);
    }

    G4ParticleDefinition* part = 
      (G4ParticleTable::GetParticleTable())->FindParticle(namePart);

    G4HadronicProcess* proc = 0;
    G4HadronicInteraction* model = 0;
    
    if ( namePart =="gamma" )
    {
       proc = new G4PhotoNuclearProcess();
    }
    
    if ( nameGen == "Bertini" )
    {
       model = new G4CascadeInterface();
    }
//    else if ( nameGen == "CHIPS" )
//    {
//       model = new G4GammaNuclearReaction();
//    }
    
    // make sure the pointers are valid
    //
    assert(proc);
    assert(model);
    
    // configure tthe model 
    // (to a specific energy range)
    //
    model->SetMinEnergy(0.9*m_p);
    model->SetMaxEnergy(1.1*m_p);
    
    // now assign model to the process
    //
    proc->RegisterMe(model);
        
    if (!proc) { 
      G4cout << "For particle: " << part->GetParticleName()
	     << " generator " << nameGen << " is unavailable"
	     << G4endl;
      exit(1);
    } 

    const G4Element* elm = material->GetElement(0); 

    G4double mass = part->GetPDGMass();
    energy = sqrt(m_p*m_p + mass*mass);

    G4cout << "energy = " << energy/GeV << " GeV" << G4endl;

    // Create a DynamicParticle
    G4DynamicParticle dParticle(part,aDirection,energy);
    
    // xsec business
    //
    G4double sigma = proc->GetElementCrossSection( &dParticle, elm );
    
    G4Navigator* nav = new G4Navigator(); // FIXME leak
    nav->SetWorldVolume(pFrame);
    G4TouchableHandle touch(nav->CreateTouchableHistory());
    
    if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
      G4cout << "G4StateManager PROBLEM! " << G4endl;

    G4Track* gTrack;
    gTrack = new G4Track(&dParticle,aTime,aPosition);
    gTrack->SetTouchableHandle(touch);

    // Step

    G4Step *step; 
    step = new G4Step(); // FIXME leak
    step->InitializeStep(gTrack);
    step->GetPreStepPoint()->SetMaterial(material);// we somehow still need it
    step->SetStepLength(theStep);
    
    // OK, we're almost done with prep work
    // final step: book histo's
    //
    std::ostringstream ossmom;
    ossmom << m_p;    
    std::string outputName = namePart + ossmom.str() + "MeV" + nameMatRaw + nameGen + ".root";
    TFile* hFile = new TFile( outputName.c_str(), "RECREATE" );
    
    TString title = namePart + "(" + ossmom.str() + "MeV) + " + nameMatRaw + " -> X ";   
//    TString xt_p    = ";Kinetic Energy [MeV]";
//    TString xyt_p   = xt_p + ";d^{2}#sigma / dE d#Omega    [#mub MeV^{-1} sr^{-1}]";
    
    G4double maxKE = dParticle.GetTotalEnergy();
    
    //
    // plots for secondary proton
    //
    TH1F* proton45  = new TH1F("p45deg", title+"+ p (45 deg)",100,0.0,maxKE);
    TH1F* proton60  = new TH1F("p60deg", title+"+ p (60 deg)",100,0.0,maxKE);
    TH1F* proton72  = new TH1F("p72deg", title+"+ p (72 deg)",100,0.0,maxKE);
    TH1F* proton84  = new TH1F("p84deg", title+"+ p (84 deg)",100,0.0,maxKE);
    TH1F* proton90  = new TH1F("p90deg", title+"+ p (90 deg)",100,0.0,maxKE);
    TH1F* proton135 = new TH1F("p135deg",title+"+ p (135 deg)",100,0.0,maxKE);
    TH1F* proKE     = new TH1F("proKE", title+"proton; Kinetic Energy (MeV)",100,0.0,maxKE);
    TH1F* proCosTh  = new TH1F("pCosTheta",  title+"p;cos #Theta",  100,-1.,1.);
    
    //
    // plots for secondary pi-
    //
    
    TH1F* pim28deg = new TH1F("pim28deg", title+"+ pi- (28.4 deg)", 100, 0., maxKE );
    pim28deg->GetXaxis()->SetTitle("momentum of secondary pi- (MeV/c)");
    TH1F* pim44deg = new TH1F("pim44deg", title+"+ pi- (44.2 deg)", 100, 0., maxKE );
    pim44deg->GetXaxis()->SetTitle("momentum of secondary pi- (MeV/c)");
    
    //
    // plots for secondary pi+
    //
    TH1F* pip28deg = new TH1F("pip28deg", title+"+ pi+ (28.4 deg)", 100, 0., maxKE );
    pip28deg->GetXaxis()->SetTitle("momentum of secondary pi+ (MeV/c)");
    TH1F* pip44deg = new TH1F("pip44deg", title+"+ pi+ (44.2 deg)", 100, 0., maxKE );
    pip44deg->GetXaxis()->SetTitle("momentum of secondary pi+ (MeV/c)");
   
    G4double weight = (sigma/microbarn)/nevt ;
    
// -> leftover from test48, for parallel processing:    if ( jobid > -1 ) histo.SetJobID(jobid);
    
    G4Timer* timer = new G4Timer();
    timer->Start();
    G4ThreeVector   mom;
    G4LorentzVector labv;
    G4VParticleChange* aChange = 0;
    //    G4VParticleChange* bChange = 0;

    if (verbose>=2) {
      G4cout << "Possible secondary status codes: " 
             << " fAlive " << fAlive 
             << " fStopButAlive " << fStopButAlive 
             << " fStopAndKill " << fStopAndKill 
             << " fKillTrackAndSecondaries " << fKillTrackAndSecondaries
             << " fSuspend " << fSuspend
             << " fPostponeToNextEvent " << fPostponeToNextEvent
             << G4endl;
    }
    
    for (G4int iter=0; iter<nevt; ++iter) 
    {

      if(verbose>=1 and iter == 1000*(iter/1000))  
        G4cout << "### " << iter << "-th event start " << G4endl;

      G4double e0 = energy-mass;

      dParticle.SetKineticEnergy(e0);

      gTrack->SetStep(step);
      gTrack->SetKineticEnergy(e0);
      gTrack->SetTouchableHandle(touch);// Set Box touchable history 
                                        // - this darn stuff is needed by CHIPS      
      aChange = proc->PostStepDoIt(*gTrack,*step);
      
      if(verbose>=1 and iter == 1000*(iter/1000)) 
        G4cout << "##### " << iter << "-th event  #####" << G4endl;
      else if (verbose>=2) 
        G4cout << "##### " << iter << "-th event  #####" << G4endl;

      // loop over secondaries and cleanup
      G4int nsecnd = aChange->GetNumberOfSecondaries();

      if(verbose>=0 and iter == 1000*(iter/1000)) 
        G4cout << "##### " << nsecnd << " nsecondaries #####" << G4endl;
      else if (verbose>=1)
        G4cout << "##### " << nsecnd << " nsecondaries #####" << G4endl;

      G4double wei = 1.;
      G4double binFactor = 1.;
      for (G4int i=0; i<nsecnd; ++i) 
      {
        G4Track*           sTrack = aChange->GetSecondary(i);
// -->        G4TrackStatus sts =  	sTrack->GetTrackStatus ();
        const G4DynamicParticle* sdp = sTrack->GetDynamicParticle();
        const G4ParticleDefinition* sdpd = sdp->GetDefinition();
        G4String sName = sdpd->GetParticleName();
	
	G4double secKE = sTrack->GetKineticEnergy();
	G4double secP  = sTrack->GetMomentum().mag();
	G4double cosTheta = sTrack->GetMomentumDirection().z();
	
	// secondary proton
	//
	if ( sName == "proton" )
	{
	   proKE->Fill(secKE,1.0);
	   proCosTh->Fill(cosTheta,1.0);
	   if (cosTheta > 0.6 && cosTheta <= 0.8) 
	   {
	      binFactor = 0.2*twopi * proton45->GetBinWidth(1);
	      wei = weight/binFactor;
	      proton45->Fill(secKE,wei);
	   } 
	   else if (cosTheta > 0.4 && cosTheta <= 0.6) 
	   {
	      binFactor = 0.2*twopi * proton60->GetBinWidth(1);
	      wei = weight/binFactor;
	      proton60->Fill(secKE,wei);
	   } 
	   else if (cosTheta > 0.2 && cosTheta <= 0.4) 
	   {
	      binFactor = 0.2*twopi * proton72->GetBinWidth(1);
	      wei = weight/binFactor;
	      proton72->Fill(secKE,wei);
	   } 
	   else if (cosTheta > 0.0 && cosTheta <= 0.2) 
	   {
	      binFactor = 0.2*twopi * proton84->GetBinWidth(1);
	      wei = weight/binFactor;
	      proton84->Fill(secKE,wei);
	   } 
	   else if (cosTheta > -0.1 && cosTheta < 0.1) 
	   {
	      binFactor = 0.2*twopi * proton90->GetBinWidth(1);
	      wei = weight/binFactor;
	      proton90->Fill(secKE,wei);
	   } 
	   else if (cosTheta > -0.6 && cosTheta <= 0.8) 
	   {
	      binFactor = 0.2*twopi * proton135->GetBinWidth(1);
	      wei = weight/binFactor;
	      proton135->Fill(secKE,wei);
	   }
	}
	//
	// secondary pi-
	//
	else if ( sName == "pi-" )
	{	   
	   if ( cosTheta > 0.7 && cosTheta <= 0.75 )
	   {
	      // pi- out at 44.2 degree (cos(44.2)=0.7169)
	      //
	      binFactor = 0.05*twopi * pim44deg->GetBinWidth(1);
	      wei = weight/binFactor;
	      pim44deg->Fill( secP, wei );
	   }
	   else if ( cosTheta > 0.85 && cosTheta <= 0.9 )
	   {
	      // pi- out at 28.4 (cos(28.4)=0.8796)
	      //
	      binFactor = 0.05*twopi * pim28deg->GetBinWidth(1);
	      wei = weight/binFactor;
	      pim28deg->Fill( secP, wei );
	   }
	}
	//
	// secondary pi+
	//
	else if ( sName == "pi+" )
	{
	   if ( cosTheta > 0.7 && cosTheta <= 0.75 )
	   {
	      // pi- out at 44.2 degree (cos(44.2)=0.7169)
	      //
	      binFactor = 0.05*twopi * pip44deg->GetBinWidth(1);
	      wei = weight/binFactor;
	      pip44deg->Fill( secP, wei );
	   }
	   else if ( cosTheta > 0.85 && cosTheta <= 0.9 )
	   {
	      // pi- out at 28.4 (cos(28.4)=0.8796)
	      //
	      binFactor = 0.05*twopi * pip28deg->GetBinWidth(1);
	      wei = weight/binFactor;
	      pip28deg->Fill( secP, wei );
	   }
	}

      } // end loop over secondaries

      // cleanup before the next event
      for (G4int i=0; i<nsecnd; ++i) 
      {
          delete aChange->GetSecondary(i);
      }
      aChange->Clear();
      
    } // event loop
     
    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;
    
    // write out histo's (per run/test)
    //    
    hFile->Write();
    hFile->Close();

  } // loop over models/beams/targets/...

  delete fin;

  delete trackManager;

  G4cout << "###### End of test #####" << G4endl;
}
