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

#ifdef G4ANALYSIS_USE
#include "DMXAnalysisManager.hh"
#endif

#include "G4SystemOfUnits.hh"
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
#include "G4ios.hh"
#include <fstream>
#include <iomanip>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXEventAction::DMXEventAction(DMXRunAction* DMXRun,
			       DMXPrimaryGeneratorAction* DMXGenerator) 
  : runAct(DMXRun),genAction(DMXGenerator)
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

  delete eventMessenger;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::BeginOfEventAction(const G4Event* evt) {


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
#ifdef G4ANALYSIS_USE
	// pass first event for histogram record:
	DMXAnalysisManager* analysis = DMXAnalysisManager::getInstance();
	analysis->HistTime(time);
#endif
    }
  }
  
    if (event_id%printModulo == 0 && P_hits > 0) {
      G4cout << "     Average light collection time: "
	     << G4BestUnit(aveTimePmtHits,"Time") << G4endl;
      G4cout << "     Number of PMT hits (photoelectron equivalent): " 
	     << P_hits << G4endl;
      //#ifdef G4ANALYSIS_USE
      // plot histograms, interactively:
      //G4bool plotevent=runAct->Getplotevent();      
      //if(plotevent) {
      //	DMXAnalysisManager* analysis = DMXAnalysisManager::getInstance();
      //	analysis->PlotHistosInter(P_hits);
      //}
      //#endif
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
void DMXEventAction::writeScintHitsToFile(void) {

  //  G4String filename="hits.out";
  G4String filename=runAct->GetsavehitsFile();
  std::ofstream hitsfile(filename, std::ios::app);
  if(!event_id) {
    std::ofstream hitsfile(filename);
    hitsfile <<"Evt     Eprim   Etot    LXe     LXeTime PMT     PMTTime Seed1           Seed2           First   Flags" 
	     << G4endl;
    hitsfile <<"#       MeV     MeV     hits    ns      hits    ns                                      hit"
	     << G4endl
	     << G4endl;
  }

  if(S_hits) {

    if(hitsfile.is_open()) {


      hitsfile << std::setiosflags(std::ios::fixed)
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
      hitsfile.close();
    }

#ifdef G4ANALYSIS_USE
    long seed1 = *seeds;
    long seed2 = *(seeds+1);    
    // pass event summary to analysis manager for booking into histos and ntple
    DMXAnalysisManager* analysis = DMXAnalysisManager::getInstance();
    analysis->analyseScintHits(event_id,energy_pri,totEnergy,S_hits,firstLXeHitTime,P_hits,aveTimePmtHits,firstParticleName,firstParticleE,gamma_ev,neutron_ev,positron_ev,electron_ev,other_ev,seed1,seed2);
#endif
  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::writePmtHitsToFile(const DMXPmtHitsCollection* hits) {

  //  G4String filename="pmt.out";
  G4String filename=runAct->GetsavepmtFile();
  std::ofstream pmtfile(filename, std::ios::app);
  G4double x; G4double y; G4double z;

  if(pmtfile.is_open()) {
    pmtfile << "Hit#    X, mm   Y, mm   Z, mm" << G4endl;       
    pmtfile << std::setiosflags(std::ios::fixed)
	    << std::setprecision(3)
	    << std::setiosflags(std::ios::left)
	    << std::setw(6);
    for (G4int i=0; i<P_hits; i++)
      {
	x = ((*hits)[i]->GetPos()).x()/mm;
	y = ((*hits)[i]->GetPos()).y()/mm;
	z = ((*hits)[i]->GetPos()).z()/mm;
	pmtfile << i << "\t"
		<< x << "\t" 
		<< y << "\t"
		<< z << G4endl;
#ifdef G4ANALYSIS_USE
  	// pass pmt hit summary to analysis manager for booking in ntples
	DMXAnalysisManager* analysis = DMXAnalysisManager::getInstance();
	analysis->analysePMTHits(event_id,i,x,y,z);
#endif
      }
    if (event_id%printModulo == 0 && P_hits > 0) 
      G4cout << "     " << P_hits << " PMT hits in " << filename << G4endl;  
    pmtfile.close();
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
