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

#include "XrayFluoEventAction.hh"

#include "XrayFluoSensorHit.hh"
#include "XrayFluoEventActionMessenger.hh"

#include "g4rw/tvordvec.h"

#ifdef G4ANALYSIS_USE
#include "XrayFluoAnalysisManager.hh"
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
XrayFluoEventAction::XrayFluoEventAction(XrayFluoAnalysisManager* aMgr)
  :drawFlag("all"),sampleCollID(-1),sensorCollID(-1),
  analysisManager(aMgr)
{
  //eventMessenger = new XraySpactrumEventActionMessenger(this);
}
#else

//This method is invoked by the G4EventManager when a G4Event object
//is sent to it

XrayFluoEventAction::XrayFluoEventAction()
:sensorCollID(-1),sampleCollID(-1),drawFlag("all"),eventMessenger(NULL),
 printModulo(1)
{
  eventMessenger = new XraySpactrumEventActionMessenger(this);
}
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoEventAction::~XrayFluoEventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoEventAction::BeginOfEventAction(const G4Event* evt)
{
  
  G4int evtNb = evt->GetEventID(); // Returns the event ID
  // if (evtNb%printModulo == 0)
  // { 
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    HepRandom::showEngineStatus();
    //  }
    
 if (sensorCollID==-1)
   {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    sensorCollID = SDman->GetCollectionID("SenCollection");
    //the pointer points to the ID number of the sensitive detector
}
 if (sampleCollID==-1)
   {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    sampleCollID = SDman->GetCollectionID("SamCollection");
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoEventAction::EndOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID(); //  Returns the event ID
  
  // extracted from hits, compute the total energy deposit (and total charged
  // track length) 

  //Returns a collection of hits for each event
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  
  //a collection of hits for the detector
 XrayFluoSensorHitsCollection* DeHC = NULL;  
  G4int n_hit = 0;
  G4double totEDet=0, totLDet=0; 
    
  if (HCE){ DeHC = (XrayFluoSensorHitsCollection*)(HCE->GetHC(sensorCollID));}

  //GetHC(G4int) is a member function of the class G4CofThisEvent
  //it returns a pointer to the SenCollection, a collection of hits
  //the argument is the number assigned by the G4SDManager and
  //it can be obtained through the G4SDManager::GetHitsCollectionID()
  //method
  //DeHC now points to SenCollection 

if (DeHC)
    {
     n_hit = DeHC->entries();
     //entries() returns the number of hits stored in DeHC
    
   for (G4int i=0;i<n_hit;i++)
	{
	  totEDet += (*DeHC)[i]->GetEdepDet(); 
	  totLDet += (*DeHC)[i]->GetTrakDet(); 
	}
     }
   
   // print this information event by event (modulo n)  	
	  
  if (evtNb%printModulo == 0) {
   
   G4cout
       << "   Detector: total energy: " << G4std::setw(7) << G4BestUnit(totEDet,"Energy")
       << "       total track length: " << G4std::setw(7) << G4BestUnit(totLDet,"Length")
       << G4endl;
       
    G4cout << "\n     " << n_hit
	   << " hits are stored in MySensorHitsCollection." << G4endl;    
  }
 HCE = evt->GetHCofThisEvent();  
//a collection of hits for the sample	     
   XrayFluoSensorHitsCollection* SaHC = NULL;
  G4int n_hit2 = 0;
  G4double totESam=0, totLSam=0; 
    
  if (HCE) SaHC = (XrayFluoSensorHitsCollection*)(HCE->GetHC(sampleCollID));
  //SaHC now points to SamCollection
  
 if (SaHC)
    {
     n_hit2 = SaHC->entries();
     for (G4int i=0;i<n_hit2;i++)
	{
	  totESam += (*SaHC)[i]->GetEdepSam(); 
	  totLSam += (*SaHC)[i]->GetTrakSam(); 
	 
           //Here we put the hit data in a an ASCII file for 
	  // later analysis  
  outFile<<G4std::setw(7) <<evtNb<<" " <<
           totESam/MeV<<" "<<
	    G4endl;
       }
// Here we call the analysis manager function for visualization

#ifdef G4ANALYSIS_USE
      analysisManager->EndOfEvent(n_hit2);
#endif 

}
   
// print this information event by event (modulo n)  	
	  
  if (evtNb%printModulo == 0) {
    
 G4cout
       << "   Sample: total energy: " << G4std::setw(7) << G4BestUnit(totESam,"Energy")
       << "       total track length: " << G4std::setw(7) << G4BestUnit(totLSam,"Length")
       << G4endl;
       
    G4cout << "\n     " << n_hit2
	   << " hits are stored in MySensorHitsCollection." << G4endl;  	  G4cout << "---> End of event: " << evtNb << G4endl;	    
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












