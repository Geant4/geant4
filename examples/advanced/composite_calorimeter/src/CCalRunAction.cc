///////////////////////////////////////////////////////////////////////////////
// File: CCalRunAction.cc
// Description: A class for providing user actions at begin and end of run
///////////////////////////////////////////////////////////////////////////////
#include "CCalRunAction.hh"
#include "CCalAnalysis.hh"

#include "globals.hh"
#include "G4Run.hh"

#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"


void CCalRunAction::BeginOfRunAction(const G4Run* aRun) {

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // A.R. Added for visualization of events.
  if ( G4VVisManager::GetConcreteInstance() ) {
    G4UImanager* UI = G4UImanager::GetUIpointer(); 
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  } 
  
  CCalAnalysis* analysis = CCalAnalysis::getInstance();
  analysis->BeginOfRun(aRun->GetRunID());
  
}


void CCalRunAction::EndOfRunAction(const G4Run* aRun) {

  G4cout << "Executing CCalRunAction" << endl;

  // A.R. Added for visualization of events.
  if (G4VVisManager::GetConcreteInstance()) {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }

  CCalAnalysis* analysis = CCalAnalysis::getInstance();
  analysis->EndOfRun(aRun->GetRunID());

}   

