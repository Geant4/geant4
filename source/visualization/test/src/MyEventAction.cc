// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyEventAction.cc,v 1.2 1999-12-15 14:54:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyEventAction.hh"

#include "MyTrackerHitsCollection.hh"
#include "MyCalorimeterHitsCollection.hh"

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

MyEventAction::MyEventAction()
{;}

MyEventAction::~MyEventAction()
{;}

void MyEventAction::BeginOfEventAction()
{;}

void MyEventAction::EndOfEventAction()
{
  const G4Event* evt = fpEventManager->GetConstCurrentEvent();

  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  G4int trackerCollID = SDman->GetCollectionID(colNam="TrackerCollection");
  G4int calorimeterCollID = SDman->GetCollectionID(colNam="CalCollection");

  G4cout << ">>> Event " << evt->GetEventID() << G4endl;

  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer)
  { n_trajectories = trajectoryContainer->entries(); }
  G4cout << "    " << n_trajectories 
       << " trajectories stored in this event." << G4endl;

  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  G4int n_hitCollection = 0;
  if(HCE)
  { n_hitCollection = HCE->GetCapacity(); }
  G4cout << "    " << n_hitCollection
       << " hitsCollections stored in this event." << G4endl;

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    //G4UImanager::GetUIpointer()->ApplyCommand("/vis~/clear/view");
    //G4UImanager::GetUIpointer()->ApplyCommand("/vis~/draw/current");
    for(G4int i=0; i<n_trajectories; i++)
    { (*(evt->GetTrajectoryContainer()))[i]->DrawTrajectory(50); }

    MyTrackerHitsCollection* THC 
      = (MyTrackerHitsCollection*)(HCE->GetHC(trackerCollID));
    if(THC) THC->DrawAllHits();
    MyCalorimeterHitsCollection* CHC
      = (MyCalorimeterHitsCollection*)(HCE->GetHC(calorimeterCollID));
    if(CHC) CHC->DrawAllHits();

    //G4UImanager::GetUIpointer()->ApplyCommand("/vis~/show/view");
  }
}



