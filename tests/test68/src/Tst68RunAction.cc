#include "Tst68RunAction.hh"
#include "globals.hh"
#include "G4Run.hh"


Tst68RunAction::~Tst68RunAction() {}


void Tst68RunAction::BeginOfRunAction( const G4Run* aRun ) {
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}


void Tst68RunAction::EndOfRunAction( const G4Run* ) {}   

