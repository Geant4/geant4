// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyEventAction.cc,v 1.3 1999-05-07 10:22:29 barrand Exp $
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

void MyEventAction::BeginOfEventAction(const G4Event*)
{;}

void MyEventAction::EndOfEventAction(const G4Event* aEvent)
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  G4int trackerCollID = SDman->GetCollectionID(colNam="TrackerCollection");
  G4int calorimeterCollID = SDman->GetCollectionID(colNam="CalCollection");

  G4cout << ">>> Event " << aEvent->GetEventID() << endl;

  G4TrajectoryContainer * trajectoryContainer = 
    aEvent->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer)
  { n_trajectories = trajectoryContainer->entries(); }
  G4cout << "    " << n_trajectories 
       << " trajectories stored in this event." << endl;

  G4HCofThisEvent * HCE = aEvent->GetHCofThisEvent();
  G4int n_hitCollection = 0;
  if(HCE)
  { n_hitCollection = HCE->GetCapacity(); }
  G4cout << "    " << n_hitCollection
       << " hitsCollections stored in this event." << endl;

  /*debug
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    //G4UImanager::GetUIpointer()->ApplyCommand("/vis~/clear/view");
    //G4UImanager::GetUIpointer()->ApplyCommand("/vis~/draw/current");
    //int mode = 50; // Draw circle of screen size 0.0001 * mode.
    int mode = 0;
    for(G4int i=0; i<n_trajectories; i++)
    { (*(aEvent->GetTrajectoryContainer()))[i]->DrawTrajectory(mode); }

    MyTrackerHitsCollection* THC 
      = (MyTrackerHitsCollection*)(HCE->GetHC(trackerCollID));
    if(THC) THC->DrawAllHits();
    MyCalorimeterHitsCollection* CHC
      = (MyCalorimeterHitsCollection*)(HCE->GetHC(calorimeterCollID));
    if(CHC) CHC->DrawAllHits();

    //G4UImanager::GetUIpointer()->ApplyCommand("/vis~/show/view");
  }
  */
}



