
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
  G4cout << ">>> Event " << evt->GetEventID() << endl;

  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  MyTrackerHitsCollection * pTHC1
   = (MyTrackerHitsCollection*)(HCE->GetHC(colID1));
  MyTrackerHitsCollection * pTHC2
   = (MyTrackerHitsCollection*)(HCE->GetHC(colID2));

  G4cout << "Event : " << evt->GetEventID() << endl;
  G4cout << pTHC1->GetSDname() << "/" << pTHC1->GetName() << endl;
  G4cout << "Number of hits          " << pTHC1->entries() << endl;
  G4cout << "  Edep of the first Hit " << (*pTHC1)[0]->GetEdep() << endl;
  G4cout << pTHC2->GetSDname() << "/" << pTHC2->GetName() << endl;
  G4cout << "Number of hits          " << pTHC2->entries() << endl;
  G4cout << "  Edep of the first Hit " << (*pTHC2)[0]->GetEdep() << endl;
}



