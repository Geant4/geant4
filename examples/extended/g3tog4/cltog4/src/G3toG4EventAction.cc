// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4EventAction.cc,v 1.1 2000-07-24 11:23:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "g4rw/tvordvec.h"
#include "G4ios.hh"
#include "G3toG4EventAction.hh"
#include "G3toG4EventActionMessenger.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G3toG4EventAction::G3toG4EventAction()
  : drawFlag("all"),eventMessenger(NULL)
{
  eventMessenger = new G3toG4EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G3toG4EventAction::~G3toG4EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G3toG4EventAction::BeginOfEventAction(const G4Event* Ev)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G3toG4EventAction::EndOfEventAction(const G4Event* Ev)
{
  const G4Event* evt = fpEventManager->GetConstCurrentEvent();

  G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  
  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer){ 
    n_trajectories = trajectoryContainer->entries(); 
  }
  G4cout << "    " << n_trajectories 
	 << " trajectories stored in this event." << G4endl;

  if(G4VVisManager::GetConcreteInstance()){
    for(G4int i=0; i<n_trajectories; i++) {
      G4Trajectory* trj = (G4Trajectory*)(*(evt->GetTrajectoryContainer()))[i];
      if (drawFlag == "all") trj->DrawTrajectory(50);
      else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
	trj->DrawTrajectory(50); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


