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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "fluoTestEventAction.hh"

#include "fluoTestSensorHit.hh"
#include "fluoTestEventActionMessenger.hh"

#include "g4rw/tvordvec.h"

#ifdef G4ANALYSIS_USE
#include "fluoTestAnalysisManager.hh"
#endif

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"


// This file is a global variable in which we store energy deposition per hit
// and other relevant information

extern G4std::ofstream outFile;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
fluoTestEventAction::fluoTestEventAction(fluoTestAnalysisManager* aMgr)
  :drawFlag("all"),sampleCollID(-1),siCollID(-1),hpgeCollID(-1),
   analysisManager(aMgr)
{
  //eventMessenger = new fluoTestEventActionMessenger(this);
}
#else


fluoTestEventAction::fluoTestEventAction()
  :siCollID(-1),sampleCollID(-1),hpgeCollID(-1),drawFlag("all"),
   eventMessenger(NULL), printModulo(1)
{
  eventMessenger = new fluoTestEventActionMessenger(this);
}
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

fluoTestEventAction::~fluoTestEventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void fluoTestEventAction::BeginOfEventAction(const G4Event* evt)
{
  
  G4int evtNb = evt->GetEventID(); // Returns the event ID
 
  G4cout << "\n---> Begin of event: " << evtNb << G4endl;
  HepRandom::showEngineStatus();
    
  if (siCollID==-1)
    {
      G4SDManager * SDman = G4SDManager::GetSDMpointer();
      siCollID = SDman->GetCollectionID("SiCollection");
     
    }
  if (sampleCollID==-1)
    {
      G4SDManager * SDman = G4SDManager::GetSDMpointer();
      sampleCollID = SDman->GetCollectionID("SamCollection");
    }
  if (hpgeCollID==-1)
    {
      G4SDManager * SDman = G4SDManager::GetSDMpointer();
      hpgeCollID = SDman->GetCollectionID("HPGeCollection");
    }
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void fluoTestEventAction::EndOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID(); //  Returns the event ID
  
  // extracted from hits, compute the total energy deposit (and total charged
  // track length) 

  //Returns a collection of hits for each event
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  
  //a collection of hits for the detector
  fluoTestSensorHitsCollection* SiHC = NULL;  
  G4int n_hit = 0;
  G4double totESi=0, totLSi=0; 
  
  if (HCE){ SiHC = (fluoTestSensorHitsCollection*)(HCE->GetHC(siCollID));}

  if (SiHC)
    {
      n_hit = SiHC->entries();
     //entries() returns the number of hits stored in DeHC
    
      for (G4int i=0;i<n_hit;i++)
	{
	  totESi += (*SiHC)[i]->GetEdepSi(); 
	  totLSi += (*SiHC)[i]->GetTrakSi(); 
	}
    }
   
  // print this information event by event (modulo n)  	
	  
  if (evtNb%printModulo == 0) {
   
    G4cout
      << "   Detector: total energy: " << G4std::setw(7) << G4BestUnit(totESi,"Energy")
       << "       total track length: " << G4std::setw(7) << G4BestUnit(totLSi,"Length")
       << G4endl;
       
    G4cout << "\n     " << n_hit
	   << " hits are stored in FluoTestSensorHitsCollection." << G4endl;    
  }
  HCE = evt->GetHCofThisEvent();  


  //a collection of hits for the sample	     
  fluoTestSensorHitsCollection* SaHC = NULL;
  G4int n_hit2 = 0;
  G4double totESam=0, totLSam=0; 
    
  if (HCE) SaHC = (fluoTestSensorHitsCollection*)(HCE->GetHC(sampleCollID));
  //SaHC now points to SamCollection
  
  if (SaHC)
    {
      n_hit2 = SaHC->entries();
      for (G4int i=0;i<n_hit2;i++)
	{
	  totESam += (*SaHC)[i]->GetEdepSam(); 
	  totLSam += (*SaHC)[i]->GetTrakSam(); 
	}
    }

   
  // print this information event by event (modulo n)  	
	  
  if (evtNb%printModulo == 0) {
    
    G4cout
      << "   Sample: total energy: " << G4std::setw(7) << G4BestUnit(totESam,"Energy")
      << "       total track length: " << G4std::setw(7) << G4BestUnit(totLSam,"Length")
      << G4endl;
       
    G4cout << "\n     " << n_hit2
	   << " hits are stored in FluoTestSensorHitsCollection." 
	   << G4endl;  	  
    G4cout << "---> End of event: " << evtNb << G4endl;	    
  }
  HCE = evt->GetHCofThisEvent();  
  //a collection of hits for the HPGeDetector     
  fluoTestSensorHitsCollection* HPGeHC = NULL;
  G4int n_hit3 = 0;
  G4double totEHPGe=0, totLHPGe=0; 
    
  if (HCE) HPGeHC = (fluoTestSensorHitsCollection*)(HCE->GetHC(hpgeCollID));
  //HPGeHC now points to SamCollection
  
  if (HPGeHC)
    {
      n_hit3 = HPGeHC->entries();
      for (G4int i=0;i<n_hit3;i++)
	{
	  totEHPGe += (*SaHC)[i]->GetEdepHPGe(); 
	  totLHPGe += (*SaHC)[i]->GetTrakHPGe(); 
	}
#ifdef G4ANALYSIS_USE
      analysisManager->EndOfEvent(n_hit3);
#endif 
    } 
	  // extract the trajectories and draw them
  
  if (G4VVisManager::GetConcreteInstance())
    {
   
      G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
      G4int n_trajectories = 0;
      if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

      for (G4int i=0; i<n_trajectories; i++) 
	{ G4Trajectory* trj = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
	if (drawFlag == "all") trj->DrawTrajectory(50);
	else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
	  trj->DrawTrajectory(50);
	else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.))
	  trj->DrawTrajectory(50);				   
	}
    }             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
