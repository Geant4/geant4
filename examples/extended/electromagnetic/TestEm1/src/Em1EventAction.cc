// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1EventAction.cc,v 1.2 1999-12-15 14:48:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em1EventAction.hh"

#include "Em1RunAction.hh"
#include "Em1EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1EventAction::Em1EventAction(Em1RunAction* run)
:Em1Run(run),drawFlag("all"),printModulo(10000),eventMessenger(NULL)
{
  eventMessenger = new Em1EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1EventAction::~Em1EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1EventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) 
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
    
 TotalEnergyDeposit = 0.;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1EventAction::EndOfEventAction(const G4Event* evt)
{
  if (drawFlag != "none") G4cout << " Energy deposit: " 
                                 << G4BestUnit(TotalEnergyDeposit,"Energy")
                                 << G4endl;

  if (G4VVisManager::GetConcreteInstance())
  {
   const G4Event* evt = fpEventManager->GetConstCurrentEvent(); 
   G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
   G4int n_trajectories = 0;
   if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();  
   for(G4int i=0; i<n_trajectories; i++) 
      { G4Trajectory* trj = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
        if (drawFlag == "all") trj->DrawTrajectory(50);
        else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                               trj->DrawTrajectory(50); 
      }
  }
  
  //save rndm status
  if (Em1Run->GetRndmFreq() == 2)
    { 
     HepRandom::saveEngineStatus("endOfEvent.rndm");   
     G4int evtNb = evt->GetEventID();
     if (evtNb%printModulo == 0)
       { 
        G4cout << "\n---> End of Event: " << evtNb << G4endl;
        HepRandom::showEngineStatus();
       }
    }     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


