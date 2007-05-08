#include "MyEventAction.hh"
#include "G4ios.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"


MyEventAction::MyEventAction() {}


MyEventAction::~MyEventAction() {}


void MyEventAction::BeginOfEventAction( const G4Event* ) { 
  //G4cout << "\n---> Begin of event: " << evt->GetEventID() << G4endl;
}


void MyEventAction::EndOfEventAction( const G4Event* evt ) {
  G4cout << " ---  MyEventAction::EndOfEventAction  ---    event = " 
	 << evt->GetEventID() << G4endl;
}

