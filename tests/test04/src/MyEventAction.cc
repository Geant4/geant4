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

#include "MyEventAction.hh"

#include "MyTrackerHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

MyEventAction::MyEventAction()
{
  colID1 = -1;
  colID2 = -1;
}

MyEventAction::~MyEventAction()
{;}

void MyEventAction::BeginOfEventAction(const G4Event*)
{
  if(colID1<0||colID2<0)
  {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    G4String colNam;
    colID1 = SDman->GetCollectionID(colNam="EvenCollection");
    colID2 = SDman->GetCollectionID(colNam="OddCollection");
  }
}

void MyEventAction::EndOfEventAction(const G4Event* evt)
{
  G4cout << ">>> Event " << evt->GetEventID() << G4endl;

  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  MyTrackerHitsCollection * pTHC1
   = (MyTrackerHitsCollection*)(HCE->GetHC(colID1));
  MyTrackerHitsCollection * pTHC2
   = (MyTrackerHitsCollection*)(HCE->GetHC(colID2));

  G4cout << "Event : " << evt->GetEventID() << G4endl;
  G4cout << pTHC1->GetSDname() << "/" << pTHC1->GetName() << G4endl;
  G4cout << "Number of hits          " << pTHC1->entries() << G4endl;
  G4cout << "  Edep of the first Hit " << (*pTHC1)[0]->GetEdep() << G4endl;
  G4cout << pTHC2->GetSDname() << "/" << pTHC2->GetName() << G4endl;
  G4cout << "Number of hits          " << pTHC2->entries() << G4endl;
  G4cout << "  Edep of the first Hit " << (*pTHC2)[0]->GetEdep() << G4endl;
}



