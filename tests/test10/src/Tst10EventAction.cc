// $Id: Tst10EventAction.cc,v 1.4 2000-01-12 15:58:28 stesting Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This is a version for maximum particle set
//	History
//        first version              09  Sept. 1998 by S.Magni
// ------------------------------------------------------------

#include "Tst10EventAction.hh"


#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "g4std/iostream"

Tst10EventAction::Tst10EventAction()
{;}

Tst10EventAction::~Tst10EventAction()
{;}

void Tst10EventAction::BeginOfEventAction(const G4Event* evt)
{
  G4cout << ">>> Start Event " << evt->GetEventID() << G4endl;
}

void Tst10EventAction::EndOfEventAction(const G4Event* evt)
{
  G4cout << ">>> End Event " << evt->GetEventID() << G4endl;

  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer)
  { n_trajectories = trajectoryContainer->entries(); }

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager) {
    for(G4int i=0; i<n_trajectories; i++)
    { (*(evt->GetTrajectoryContainer()))[i]->DrawTrajectory(50); }
	}		 
}



