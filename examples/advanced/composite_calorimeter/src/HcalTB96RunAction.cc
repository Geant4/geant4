#include "HcalTB96RunAction.hh"
#include "HcalTB96Analysis.hh"

#include "globals.hh"
#include "G4Run.hh"

void HcalTB96RunAction::StartOfRunAction(const G4Run* aRun) {

  HcalTB96Analysis* analysis = HcalTB96Analysis::getInstance();
  analysis->BeginOfRun(aRun->GetRunID());
}

void HcalTB96RunAction::EndOfRunAction(const G4Run* aRun) {

  cout << "Executing HcalTB96RunAction" << endl;
  HcalTB96Analysis* analysis = HcalTB96Analysis::getInstance();
  analysis->EndOfRun(aRun->GetRunID());

}   
