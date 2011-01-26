
#include "RunAction.hh"
#include "AnalysisManager.hh"
#include "G4Run.hh"


void RunAction::BeginOfRunAction(const G4Run* run) {

  std::cout << "INFORMATION: Run No " << run -> GetRunID() 
            << " starts." << std::endl;
}


void RunAction::EndOfRunAction(const G4Run* run) {

  AnalysisManager::Instance() -> PrintResults();
  AnalysisManager::Instance() -> Destroy();

  std::cout << "INFORMATION: Run No " << run -> GetRunID() 
            << " ends (Number of events = " 
            << run -> GetNumberOfEvent() << ")." << std::endl;
}
