///////////////////////////////////////////////////////////////////////////////
// File: CCalRunAction.cc
// Description: A class for providing user actions at begin and end of run
///////////////////////////////////////////////////////////////////////////////
#include "CCalRunAction.hh"
#include "CCalAnalysis.hh"

#include "globals.hh"
#include "G4Run.hh"

void CCalRunAction::BeginOfRunAction(const G4Run* aRun) {

  CCalAnalysis* analysis = CCalAnalysis::getInstance();
  analysis->BeginOfRun(aRun->GetRunID());
}

void CCalRunAction::EndOfRunAction(const G4Run* aRun) {

  cout << "Executing CCalRunAction" << endl;
  CCalAnalysis* analysis = CCalAnalysis::getInstance();
  analysis->EndOfRun(aRun->GetRunID());

}   
