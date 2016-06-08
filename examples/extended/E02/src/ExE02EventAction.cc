
#include "ExE02EventAction.hh"

#include "ExE02TrackerHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

ExE02EventAction::ExE02EventAction()
{
  colID1 = -1;
  colID2 = -1;
}

ExE02EventAction::~ExE02EventAction()
{;}

void ExE02EventAction::BeginOfEventAction(const G4Event* )
{
  if(colID1<0||colID2<0)
  {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    G4String colNam;
    colID1 = SDman->GetCollectionID(colNam="EvenCollection");
    colID2 = SDman->GetCollectionID(colNam="OddCollection");
  }
}

void ExE02EventAction::EndOfEventAction(const G4Event* evt)
{
  G4cout << ">>> Event " << evt->GetEventID() << endl;

  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  ExE02TrackerHitsCollection * pTHC1
   = (ExE02TrackerHitsCollection*)(HCE->GetHC(colID1));
  ExE02TrackerHitsCollection * pTHC2
   = (ExE02TrackerHitsCollection*)(HCE->GetHC(colID2));

  G4cout << "Event : " << evt->GetEventID() << endl;
  G4cout << pTHC1->GetSDname() << "/" << pTHC1->GetName() << endl;
  G4cout << "Number of hits          " << pTHC1->entries() << endl;
  G4cout << "  Edep of the first Hit " << (*pTHC1)[0]->GetEdep() << endl;
  G4cout << pTHC2->GetSDname() << "/" << pTHC2->GetName() << endl;
  G4cout << "Number of hits          " << pTHC2->entries() << endl;
  G4cout << "  Edep of the first Hit " << (*pTHC2)[0]->GetEdep() << endl;
}



