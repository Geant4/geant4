#include "StatAccepTestRunAction.hh"
#include "globals.hh"
#include "G4Run.hh"

#include "StatAccepTestAnalysis.hh"


StatAccepTestRunAction::~StatAccepTestRunAction() { 
  // Close the AIDA file at the end of the job (not at the end of the run!).
  StatAccepTestAnalysis::getInstance()->close();  
}


void StatAccepTestRunAction::BeginOfRunAction( const G4Run* aRun ) {
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  // Reset the histograms/ntuples of the eventual previous run.
  StatAccepTestAnalysis::getInstance()->init();
}


void StatAccepTestRunAction::EndOfRunAction( const G4Run* ) {
  // Commit the histograms/ntuples.
  StatAccepTestAnalysis::getInstance()->finish();
}   

