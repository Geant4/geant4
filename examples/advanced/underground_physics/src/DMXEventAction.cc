//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// EventAction program
// --------------------------------------------------------------

#include "DMXEventAction.hh"

// note DMXPmtHit.hh and DMXScintHit.hh are included in DMXEventAction.hh

#include "DMXEventActionMessenger.hh"

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
#include "g4std/fstream"
#include "g4std/iomanip"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXEventAction::DMXEventAction()  {

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

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXEventAction::~DMXEventAction() {

  delete eventMessenger;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::BeginOfEventAction(const G4Event* evt) {

  event_id = evt->GetEventID();
 
  // print this information event by event (modulo n)  	
  if (event_id%printModulo == 0)
    G4cout << "\n---> Begin of event: " << event_id << G4endl;

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
  hitTime           = 0.;
  aveTimeScintHits  = 0.;
  particleEnergy    = 0.;
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
	if (event_id%printModulo == 0)
	  G4cout << "     First hit: " << firstParticleName << G4endl;
      }
      hitEnergy         = (*SHC)[i]->GetEdep();
      totEnergy        += hitEnergy;
      hitTime           = (*SHC)[i]->GetTime();
      aveTimeScintHits += hitTime/(double)S_hits;
      particleName      = (*SHC)[i]->GetParticle();
      particleEnergy    = (*SHC)[i]->GetParticleEnergy();

      if(particleName == "gamma") {
	gamma_ev = true;
	start_gamma = true;
	start_neutron = false;
	// globalTime = hitTime;
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
    for (G4int i=0; i<P_hits; i++)
      aveTimePmtHits += ((*PHC)[i]->GetTime())/(double)P_hits;
    if (event_id%printModulo == 0)
      G4cout << "     Average light collection time: "
	     << aveTimePmtHits / ns << " ns" << G4endl;

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
    G4cout << "---> End of event: " << event_id << G4endl;	


}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::writeScintHitsToFile(void) {

  G4String filename="analysis/hits.out";
  G4std::ofstream hitsfile(filename, G4std::ios::app);
  if(!event_id) {
    G4std::ofstream hitsfile(filename);
    hitsfile <<"Evt     Etot    LXe     LXeTime PMT     PmtTime First   Flags" 
	     << G4endl;
    hitsfile <<"#       MeV     hits    ns      hits    ns      hit" << G4endl
	     << G4endl;
  }


  if(hitsfile.is_open()) {

    hitsfile << G4std::setiosflags(G4std::ios::fixed)
             << G4std::setprecision(4)
             << G4std::setiosflags(G4std::ios::left)
             << G4std::setw(6);

    hitsfile << event_id << "\t"
             << totEnergy/MeV << "\t"
             << S_hits  << "\t"
             << aveTimeScintHits/nanosecond << "\t"
             << P_hits << "\t"
             << aveTimePmtHits/nanosecond << "\t"
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


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::writePmtHitsToFile(const DMXPmtHitsCollection* hits) {

  G4String filename="analysis/pmt.out";
  G4std::ofstream pmtfile(filename);

  if(pmtfile.is_open()) {
    pmtfile << "Hit#    X, mm   Y, mm   Z, mm" << G4endl;       
    pmtfile << G4std::setiosflags(G4std::ios::fixed)
            << G4std::setprecision(3)
            << G4std::setiosflags(G4std::ios::left)
            << G4std::setw(6);
    for (G4int i=0; i<P_hits; i++)
      pmtfile << i << "\t"
              << ((*hits)[i]->GetPos()).x()/mm << "\t" 
              << ((*hits)[i]->GetPos()).y()/mm << "\t"
              << ((*hits)[i]->GetPos()).z()/mm << G4endl;
    pmtfile.close();
    if (event_id%printModulo == 0)
      G4cout << "     PMT hits in " << filename << G4endl;  
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
