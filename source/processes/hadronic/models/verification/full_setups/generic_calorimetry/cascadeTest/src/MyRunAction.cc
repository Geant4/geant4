#include "MyRunAction.hh"
#include "globals.hh"
#include "G4Run.hh"

#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

#include "MyAnalysis.hh"


void MyRunAction::BeginOfRunAction(const G4Run* aRun) {

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // Reset the histograms/ntuples of the eventual previous run.
  MyAnalysis::getInstance()->init();

  // For the visualization of the event.
  if ( G4VVisManager::GetConcreteInstance() ) {
    G4UImanager* UI = G4UImanager::GetUIpointer(); 
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  } 
  
}


void MyRunAction::EndOfRunAction(const G4Run* aRun) {

  // Commit the histograms/ntuples.
  MyAnalysis::getInstance()->finish();

  // For the visualization of the event.
  if (G4VVisManager::GetConcreteInstance()) {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }

}   

