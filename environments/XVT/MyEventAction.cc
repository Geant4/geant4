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
// $Id: MyEventAction.cc,v 1.4 2001-07-11 09:56:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyEventAction.hh"

#include "MyCalorimeterHit.hh"
#include "MyCalorimeterHitsCollection.hh"
#include "g4rw/tvordvec.h"

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

#include "MyEventActionMessenger.cc"

MyEventAction::MyEventAction()
:drawFlag(false)
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;

  calorimeterCollID = SDman->GetCollectionID(colNam="CalCollection");
  //////////bulkCollID = SDman->GetCollectionID(colNam="bulkCollection");

  new MyEventActionMessenger(this);
}

MyEventAction::~MyEventAction()
{;}

void MyEventAction::BeginOfEventAction()
{
  if(drawFlag)
  {
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager)
    {
      G4UImanager::GetUIpointer()->ApplyCommand("/vis~/Draw/current");
    }
  }
}

void MyEventAction::EndOfEventAction()
{
  const G4Event* evt = fpEventManager->GetConstCurrentEvent();

  G4cout << ">>> Event " << evt->GetEventID() << G4endl;

  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  MyCalorimeterHitsCollection* CHC = NULL;
  if(HCE)
    CHC = (MyCalorimeterHitsCollection*)(HCE->GetHC(calorimeterCollID));

  if(CHC)
  {
    G4RWTValOrderedVector<MyCalorimeterHit> theCollection
      = CHC->GetVector();
    int n_hit = theCollection.entries();
    G4cout << "     " << n_hit
         << " hits are stored in MyCalorimeterHitsCollection." << G4endl;
    G4double totE = 0;
    for(int i=0;i<n_hit;i++)
    { totE += theCollection[i].GetEdep(); }
    G4cout << "     Total energy deposition in calorimeter tubes : "
         << totE / GeV << " (GeV)" << G4endl;
  }

  if(drawFlag)
  {
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager)
    {
      if(CHC) CHC->DrawAllHits();
      G4UImanager::GetUIpointer()->ApplyCommand("/vis~/show/view");
    }
  }
}



