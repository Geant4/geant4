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
// $Id: MyEventAction.cc,v 1.4 2001-08-24 20:25:59 johna Exp $
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
#include "G4VisAttributes.hh"
#include "G4Scale.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

//#define DRAWTRAJHIT

MyEventAction::MyEventAction()
{;}

MyEventAction::~MyEventAction()
{;}

void MyEventAction::BeginOfEventAction(const G4Event* anEvent)
{;}

void MyEventAction::EndOfEventAction(const G4Event* anEvent)
{
  static int coutCount = 0;
  if (coutCount < 10) {
    coutCount++;
    G4cout << "MySteppingAction::EndOfEventActionAction called." << G4endl;
  }

  const G4Event* evt = anEvent;

  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  G4int trackerCollID = SDman->GetCollectionID(colNam="TrackerCollection");
  G4int calorimeterCollID = SDman->GetCollectionID(colNam="CalCollection");

  if (coutCount < 10) {
    G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  }

#ifdef DRAWTRAJHIT

  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer)
  { n_trajectories = trajectoryContainer->entries(); }
  if (coutCount < 10) {
    G4cout << "    " << n_trajectories 
	   << " trajectories stored in this event." << G4endl;
  }

  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  G4int n_hitCollection = 0;
  if(HCE)
  { n_hitCollection = HCE->GetCapacity(); }
  if (coutCount < 10) {
    G4cout << "    " << n_hitCollection
	   << " hitsCollections stored in this event." << G4endl;
  }

#endif

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {

#ifdef DRAWTRAJHIT

    for(G4int i=0; i<n_trajectories; i++)
    { (*(evt->GetTrajectoryContainer()))[i]->DrawTrajectory(50); }

    MyTrackerHitsCollection* THC 
      = (MyTrackerHitsCollection*)(HCE->GetHC(trackerCollID));
    if(THC) THC->DrawAllHits();
    MyCalorimeterHitsCollection* CHC
      = (MyCalorimeterHitsCollection*)(HCE->GetHC(calorimeterCollID));
    if(CHC) CHC->DrawAllHits();

#endif

    G4Scale scale(1. * m, "Test Scale");
    G4VisAttributes va (G4Colour(1.,0.,0.));
    scale.SetVisAttributes(&va);
    pVVisManager->Draw(scale);

  }
}



