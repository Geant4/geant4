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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// History/Additions:
// 16 Jan 2002  Added analysis
//
//
// EventAction program
// --------------------------------------------------------------

#include "DMXEventAction.hh"

// pass parameters for messengers:
#include "DMXRunAction.hh"
#include "DMXPrimaryGeneratorAction.hh"

// note DMXPmtHit.hh and DMXScintHit.hh are included in DMXEventAction.hh

#include "DMXEventActionMessenger.hh"
#include "DMXAnalysisManager.hh"


#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"
#include <fstream>
#include <iomanip>

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXEventAction::DMXEventAction() 
  : runAct(0),genAction(0),hitsfile(0),pmtfile(0)
{

  // create messenger
  eventMessenger = new DMXEventActionMessenger(this);

  // defaults for messenger
  drawColsFlag = "standard";
  drawTrksFlag = "all";
  drawHitsFlag = 1;
  savePmtFlag  = 0;
  saveHitsFlag = 1;
  
  printModulo = 1;

  // hits collections
  scintillatorCollID = -1;
  pmtCollID = -1;

  energy_pri=0;
  seeds=NULL;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXEventAction::~DMXEventAction() {

  if (hitsfile)
    {
      hitsfile->close();
      delete hitsfile;
    }
  if (pmtfile)
    {
      pmtfile->close();
      delete pmtfile;
    }
  delete eventMessenger;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::BeginOfEventAction(const G4Event* evt) 
{

  //thread-local run action
  if (!runAct) 
    runAct = 
      dynamic_cast<const DMXRunAction*>
      (G4RunManager::GetRunManager()->GetUserRunAction());
  
  if (!genAction)
    genAction = dynamic_cast<const DMXPrimaryGeneratorAction*>
      (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());


  // grab seeds
  seeds = genAction->GetEventSeeds();

  // grab energy of primary
  energy_pri = genAction->GetEnergyPrimary();

  event_id = evt->GetEventID();
 
  // print this information event by event (modulo n)  	
  if (event_id%printModulo == 0)
    {
      G4cout << "\n---> Begin of event: " << event_id << G4endl;
      G4cout << "\n     Primary Energy: " << G4BestUnit(energy_pri,"Energy") 
	     << G4endl;
      //      HepRandom::showEngineStatus(); 
    }


  // get ID for scintillator hits collection
  if (scintillatorCollID==-1) {
    G4SDManager *SDman = G4SDManager::GetSDMpointer();
    scintillatorCollID = SDman->GetCollectionID("scintillatorCollection");
  }

  // get ID for pmt hits collection
  if (pmtCollID==-1) {
    G4SDManager *SDman = G4SDManager::GetSDMpointer();
    pmtCollID = SDman->GetCollectionID("pmtCollection");
  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::EndOfEventAction(const G4Event* evt) {

  // check that both hits collections have been defined
  if(scintillatorCollID<0||pmtCollID<0) return;

  G4AnalysisManager* man = G4AnalysisManager::Instance();

  // address hits collections
  DMXScintHitsCollection* SHC = NULL;
  DMXPmtHitsCollection*   PHC = NULL;
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(HCE) {
    SHC = (DMXScintHitsCollection*)(HCE->GetHC(scintillatorCollID));
    PHC = (DMXPmtHitsCollection*)(HCE->GetHC(pmtCollID));
  }

  // event summary
  totEnergy         = 0.;
  totEnergyGammas   = 0.;
  totEnergyNeutrons = 0.;
  firstParticleE    = 0.;
  particleEnergy    = 0.;
  firstLXeHitTime   = 0.;
  aveTimePmtHits    = 0.;

  firstParticleName = "";
  particleName      = "";


  // particle flags
  gamma_ev          = false;
  neutron_ev        = false;
  positron_ev       = false;
  electron_ev       = false;
  proton_ev         = false;
  other_ev          = false;
  start_gamma       = false;
  start_neutron     = false;


  // scintillator hits
  if(SHC) {
    S_hits = SHC->entries();
    
    for (G4int i=0; i<S_hits; i++) {
      if(i==0) {
	firstParticleName = (*SHC)[0]->GetParticle();
	firstLXeHitTime   = (*SHC)[0]->GetTime();
	firstParticleE = (*SHC)[0]->GetParticleEnergy();
	if (event_id%printModulo == 0 && S_hits > 0) {
	  G4cout << "     First hit in LXe: " << firstParticleName << G4endl;
	  G4cout << "     Number of hits in LXe: " << S_hits << G4endl; 
	}
      }
      hitEnergy         = (*SHC)[i]->GetEdep();
      totEnergy        += hitEnergy;
      
      particleName      = (*SHC)[i]->GetParticle();
      particleEnergy    = (*SHC)[i]->GetParticleEnergy();

      if(particleName == "gamma") {
	gamma_ev = true;
	start_gamma = true;
	start_neutron = false;
      }
      else if(particleName == "neutron") 
	neutron_ev = true;
      else if(particleName == "e+") 
	positron_ev = true;
      else if(particleName == "e-") 
	electron_ev = true;
      else if(particleName == "proton") 
	proton_ev = true;
      else {
	other_ev = true;
	start_gamma = false;
	start_neutron = true;
      }

      if(start_gamma && !start_neutron) 
	totEnergyGammas += hitEnergy;
      if(start_neutron && !start_gamma) 
	totEnergyNeutrons += hitEnergy;
    }
    
    if (event_id%printModulo == 0)
      G4cout << "     Total energy in LXe: " 
	     << G4BestUnit(totEnergy,"Energy") << G4endl;  
    
  }
  
  
  // PMT hits
  if(PHC) {
    P_hits = PHC->entries();
    
    // average time of PMT hits
    for (G4int i=0; i<P_hits; i++) {
      G4double time = ( (*PHC)[i]->GetTime() - firstLXeHitTime );
      aveTimePmtHits += time / (G4double)P_hits;
      ////      if (event_id == 0) {
      if(P_hits >= 0) {
	man->FillH1(7,time);
      }
    }
  
    if (event_id%printModulo == 0 && P_hits > 0) {
      G4cout << "     Average light collection time: "
	     << G4BestUnit(aveTimePmtHits,"Time") << G4endl;
      G4cout << "     Number of PMT hits (photoelectron equivalent): " 
	     << P_hits << G4endl;     
    }
    // write out (x,y,z) of PMT hits
    if (savePmtFlag)
      writePmtHitsToFile(PHC);
  }
  

  // write out event summary
  if(saveHitsFlag) 
    writeScintHitsToFile();
  
  // draw trajectories
  if(drawColsFlag=="standard" && drawTrksFlag!="none")
    drawTracks(evt);

  // hits in PMT
  if(drawHitsFlag && PHC) 
    PHC->DrawAllHits();

  // print this event by event (modulo n)  	
  if (event_id%printModulo == 0) 
    G4cout << "\n---> End of event: " << event_id << G4endl << G4endl;	

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::writeScintHitsToFile() 
{

  G4String filename="hits.out";
  if (runAct)
    filename=runAct->GetsavehitsFile();

 
  
  //First time it is inkoved
  if (!hitsfile)
    {
      //check that we are in a worker: returns -1 in a master and -2 in sequential
      //one file per thread is produced ending with ".N", with N= thread number
      if (G4Threading::G4GetThreadId() >= 0)
	{
	  std::stringstream sss;
	  sss << filename.c_str() << "." << G4Threading::G4GetThreadId();	 
	  filename = sss.str();
	  //G4cout << "Filename is: " << filename << G4endl;
	}
      
      hitsfile = new std::ofstream;
      hitsfile->open(filename);
      (*hitsfile) <<"Evt     Eprim   Etot    LXe     LXeTime PMT     PMTTime Seed1           Seed2           First   Flags" 
	       << G4endl;
      (*hitsfile) <<"#       MeV     MeV     hits    ns      hits    ns                                      hit"
	       << G4endl
	       << G4endl;
    }

  if(S_hits) {

    if(hitsfile->is_open()) {


      (*hitsfile) << std::setiosflags(std::ios::fixed)
		  << std::setprecision(4)
		  << std::setiosflags(std::ios::left)
		  << std::setw(6)
		  << event_id << "\t"
		  << energy_pri/MeV << "\t" 
		  << totEnergy/MeV << "\t"
		  << S_hits  << "\t"
		  << std::setiosflags(std::ios::scientific) 
		  << std::setprecision(2)
		  << firstLXeHitTime/nanosecond << "\t"
		  << P_hits << "\t"
		  << std::setiosflags(std::ios::fixed) 
		  << std::setprecision(4)
		  << aveTimePmtHits/nanosecond << "\t"
		  << *seeds     << "\t"
		  << *(seeds+1) << "\t"
		  << firstParticleName << "\t"
		  << (gamma_ev    ? "gamma " : "") 
		  << (neutron_ev  ? "neutron " : "") 
		  << (positron_ev ? "positron " : "") 
		  << (electron_ev ? "electron " : "") 
		  << (other_ev    ? "other " : "") 
		  << G4endl;

      if (event_id%printModulo == 0)
	G4cout << "     Event summary in file " << filename << G4endl;  
    }

    G4AnalysisManager* man = G4AnalysisManager::Instance();
    G4int firstparticleIndex = 0;
    if(firstParticleName == "gamma") firstparticleIndex = 1;
    else if (firstParticleName == "neutron") firstparticleIndex = 2;
    else if(firstParticleName == "electron") firstparticleIndex = 3;
    else if(firstParticleName == "positron") firstparticleIndex = 4;
    else{
      firstparticleIndex = 5;
      man->FillH1(3,totEnergy);
    }

    man->FillH1(4,P_hits,10); //weight
    man->FillH1(5,P_hits);

    man->FillH1(1,energy_pri/keV);
    man->FillH1(2,totEnergy/keV);
    man->FillH1(6,aveTimePmtHits/ns);

    long seed1 = *seeds;
    long seed2 = *(seeds+1);    
    
    //Fill ntuple #2
    man->FillNtupleDColumn(2,0,event_id);
    man->FillNtupleDColumn(2,1,energy_pri/keV);
    man->FillNtupleDColumn(2,2,totEnergy);
    man->FillNtupleDColumn(2,3,S_hits);
    man->FillNtupleDColumn(2,4,firstLXeHitTime);
    man->FillNtupleDColumn(2,5,P_hits);
    man->FillNtupleDColumn(2,6,aveTimePmtHits);
    man->FillNtupleDColumn(2,7,firstparticleIndex);
    man->FillNtupleDColumn(2,8,firstParticleE);
    man->FillNtupleDColumn(2,9,gamma_ev);
    man->FillNtupleDColumn(2,10,neutron_ev);
    man->FillNtupleDColumn(2,11,positron_ev);
    man->FillNtupleDColumn(2,12,electron_ev);
    man->FillNtupleDColumn(2,13,other_ev);
    man->FillNtupleDColumn(2,14,seed1);
    man->FillNtupleDColumn(2,15,seed2);
    man->AddNtupleRow(2);
    
  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::writePmtHitsToFile(const DMXPmtHitsCollection* hits) 
{ 
  G4String filename="pmt.out";
  if (runAct)
    filename=runAct->GetsavepmtFile();
  

  //first time it is invoked: create it
  if (!pmtfile)
    {
      //check that we are in a worker: returns -1 in a master and -2 in sequential
      //one file per thread is produced ending with ".N", with N= thread number
      if (G4Threading::G4GetThreadId() >= 0)
	{
	  std::stringstream sss;
	  sss << filename.c_str() << "." << G4Threading::G4GetThreadId();	 
	  filename = sss.str();
	  //G4cout << "Filename is: " << filename << G4endl;
	}
      pmtfile = new std::ofstream;
      pmtfile->open(filename);
    }


  G4double x; G4double y; G4double z;
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  if(pmtfile->is_open()) {
    (*pmtfile) << "Hit#    X, mm   Y, mm   Z, mm" << G4endl;       
    (*pmtfile) << std::setiosflags(std::ios::fixed)
	       << std::setprecision(3)
	       << std::setiosflags(std::ios::left)
	       << std::setw(6);
    for (G4int i=0; i<P_hits; i++)
      {
	x = ((*hits)[i]->GetPos()).x()/mm;
	y = ((*hits)[i]->GetPos()).y()/mm;
	z = ((*hits)[i]->GetPos()).z()/mm;
	(*pmtfile) << i << "\t"
		   << x << "\t" 
		   << y << "\t"
		   << z << G4endl;
	
	man->FillH2(1,x/mm, y/mm);  // fill(x,y)
	if (event_id == 0 ) {
	  man->FillH2(2,x,y);
	}

	//Fill ntuple #3
	man->FillNtupleDColumn(3,0,event_id);
	man->FillNtupleDColumn(3,1,(G4double) i);
	man->FillNtupleDColumn(3,2,x);
	man->FillNtupleDColumn(3,3,y);
	man->FillNtupleDColumn(3,4,z);
	man->AddNtupleRow(3);

      }
    if (event_id%printModulo == 0 && P_hits > 0) 
      G4cout << "     " << P_hits << " PMT hits in " << filename << G4endl;  
  }
  
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::drawTracks(const G4Event* evt) {

  if(G4VVisManager::GetConcreteInstance()) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");    
    G4TrajectoryContainer* trajContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;

    if(trajContainer) n_trajectories = trajContainer->entries();
    for (G4int i=0; i<n_trajectories; i++) {
      G4Trajectory* trj = (G4Trajectory*)(*trajContainer)[i];
      if (drawTrksFlag == "all") 
	trj->DrawTrajectory();
      else if ((drawTrksFlag == "charged") && (trj->GetCharge() != 0.))
	trj->DrawTrajectory();
      else if ((drawTrksFlag == "noscint") 
	       && (trj->GetParticleName() != "opticalphoton"))
	trj->DrawTrajectory();
    }
    
    // G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");    
  } 

}
