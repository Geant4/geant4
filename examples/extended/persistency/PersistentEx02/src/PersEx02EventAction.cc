// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersEx02EventAction.cc,v 1.4 1999/12/15 14:49:17 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//

#include "PersEx02EventAction.hh"

#include "PersEx02TrackerHit.hh"
#include "PersEx02TrackerHitsCollection.hh"

#include "G4PersistencyManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4PHCofThisEvent.hh"
#include "G4PersistentHitMan.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "HepODBMS/odbms/HepODBMS.h"

PersEx02EventAction::PersEx02EventAction()
{
  colID1 = -1;
  colID2 = -1;
}

PersEx02EventAction::~PersEx02EventAction()
{;}

void PersEx02EventAction::BeginOfEventAction(const G4Event* )
{
  if(colID1<0||colID2<0)
  {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    G4String colNam;
    colID1 = SDman->GetCollectionID(colNam="EvenCollection");
    colID2 = SDman->GetCollectionID(colNam="OddCollection");
  }
}

void PersEx02EventAction::EndOfEventAction(const G4Event* evt)
{
  G4cout << ">>> Event " << evt->GetEventID() << G4endl;

  G4PHCofThisEvent* pHCE = G4PersistentHitMan::GetPersistentHitMan()
                                                ->GetCurrentPHCofThisEvent();
  HepRef(PersEx02TrackerHitsCollection) pTHC1
   = (HepRef(PersEx02TrackerHitsCollection))(pHCE->GetHC(colID1));
  HepRef(PersEx02TrackerHitsCollection) pTHC2
   = (HepRef(PersEx02TrackerHitsCollection))(pHCE->GetHC(colID2));

  G4cout << "Event : " << evt->GetEventID() << G4endl
         << pTHC1->GetSDname() << "/" << pTHC1->GetName() << "    "
         << "Number of hits       " << pTHC1->size() << G4endl;
  if( pTHC1->size() > 0 )
  {
    HepRef(PersEx02TrackerHit) aHit = pTHC1->element(0);
    G4cout << "  Edep of the first Hit " << aHit->GetEdep() << G4endl;
  }
  G4cout << pTHC2->GetSDname() << "/" << pTHC2->GetName() << "     "
         << "Number of hits       " << pTHC2->size() << G4endl;
  if( pTHC2->size() > 0 )
  {
    HepRef(PersEx02TrackerHit) aHit = pTHC2->element(0);
    assert(aHit!=NULL);
    G4cout << "  Edep of the first Hit " << aHit->GetEdep() << G4endl;
  }

  G4PersistencyManager* persM
    = G4PersistencyManager::GetPersistencyManager();
  if( persM && (evt->GetEventID() % 10 == 0))
  {
    // commit and resume the sustained event DB transaction for every 10 events
    persM->Commit(kEventDB, true);
    persM->StartTransaction(kEventDB, kUpdate, true);
  }

}

