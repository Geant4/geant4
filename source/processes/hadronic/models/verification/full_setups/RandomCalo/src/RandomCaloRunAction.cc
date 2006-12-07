#include "RandomCaloRunAction.hh"
#include "globals.hh"
#include "G4Run.hh"



RandomCaloRunAction::~RandomCaloRunAction() { 
}


void RandomCaloRunAction::BeginOfRunAction( const G4Run* aRun ) {
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}


void RandomCaloRunAction::EndOfRunAction( const G4Run* ) {
}   

