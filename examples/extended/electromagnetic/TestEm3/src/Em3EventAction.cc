// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3EventAction.cc,v 1.1 1999-10-11 16:55:52 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em3EventAction.hh"

#include "Em3RunAction.hh"
#include "Em3DetectorConstruction.hh"
#include "Em3CalorHit.hh"
#include "Em3EventActionMessenger.hh"

#include <rw/tvordvec.h>

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
#include "Randomize.hh"

#include "CLHEP/Hist/HBookFile.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3EventAction::Em3EventAction(Em3RunAction* run, Em3DetectorConstruction* det)
:Em3Run(run),Detector(det),calorimeterCollID(-1),drawFlag("all"),
 eventMessenger(NULL),printModulo(10000)
{
  eventMessenger = new Em3EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3EventAction::~Em3EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3EventAction::BeginOfEventAction(const G4Event* evt)
{   
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) 
    G4cout << "\n---> Begin Of Event: " << evtNb << endl;
    
 if (calorimeterCollID==-1)
  {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    calorimeterCollID = SDman->GetCollectionID("CalCollection");
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3EventAction::EndOfEventAction(const G4Event* evt)
{ 
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  Em3CalorHitsCollection* CHC = NULL;
  G4int NbHits=0;
  G4int NbOfAbsor=Detector->GetNbOfAbsor();
  G4double totEAbs, totLAbs;
  HepTuple* ntuple=Em3Run->GetnTuple();
  char str1[6], str2[6];
  strcpy(str1,"EAbs");strcpy(str2,"LAbs");
    
  if (HCE) CHC = (Em3CalorHitsCollection*)(HCE->GetHC(calorimeterCollID));

  if (CHC)
    {
     NbHits = CHC->entries();
     for (G4int k=0; k<NbOfAbsor; k++)
        {totEAbs=totLAbs=0.;
         for (G4int j=0;j<NbHits;j++)
            {
	     totEAbs += (*CHC)[j]->GetEdepAbs(k); 
             totLAbs += (*CHC)[j]->GetTrakAbs(k);     
            }
         Em3Run->fillPerEvent(k,totEAbs,totLAbs);
      
         //fill ntuple
         //	 
	 str1[4] = str2[4] = (char)((int)('0') + k);
         ntuple->column(str1,totEAbs);
         ntuple->column(str2,totLAbs);
       }
      ntuple->dumpData();       	  
    }
    
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
        }
  } 
  
  //save rndm status
  if (Em3Run->GetRndmFreq() == 2)
    { 
     HepRandom::saveEngineStatus("endOfEvent.rndm");   
     G4int evtNb = evt->GetEventID();
     if (evtNb%printModulo == 0)
       { 
        G4cout << "\n---> End of Event: " << evtNb << endl;
        HepRandom::showEngineStatus();
       }
    }     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


