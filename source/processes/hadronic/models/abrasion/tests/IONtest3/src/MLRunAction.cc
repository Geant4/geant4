////////////////////////////////////////////////////////////////////////////////
//
#include "MLRunAction.hh"
#include "MLAnalysisManager.hh"

#include <stdlib.h>
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
////////////////////////////////////////////////////////////////////////////////
//
MLRunAction::MLRunAction ()
{}
////////////////////////////////////////////////////////////////////////////////
//
MLRunAction::~MLRunAction ()
{}
////////////////////////////////////////////////////////////////////////////////
//
void MLRunAction::BeginOfRunAction (const G4Run* )
{
// Open the file for the tracks of this run
// sprintf(name,"Tracks_%d.dat", aRun->GetRunID());
//outFile.open(name);

  // Prepare the visualization
  if (G4VVisManager::GetConcreteInstance()) {
    G4UImanager* UI = G4UImanager::GetUIpointer();
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  }

  // If analysis is used reset the histograms
  MLAnalysisManager* analysisManager = MLAnalysisManager::getInstance();
  analysisManager->BeginOfRunAction();
}
////////////////////////////////////////////////////////////////////////////////
//
void MLRunAction::EndOfRunAction (const G4Run* aRun)
{
// Run ended, update the visualization
  if (G4VVisManager::GetConcreteInstance()) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }

  // Close the file with the hits information
  //  outFile.close();

  // If analysis is used, print out the histograms
  MLAnalysisManager* analysisManager = MLAnalysisManager::getInstance();
  analysisManager->EndOfRunAction(G4double(aRun-> GetNumberOfEvent()));
}
////////////////////////////////////////////////////////////////////////////////
