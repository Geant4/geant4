#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "MCTruthTrackingAction.hh"
#include "MCTruthEventAction.hh"

#include "MCTruthManager.hh"

int main()
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  runManager->SetUserInitialization(new DetectorConstruction);
  runManager->SetUserInitialization(new PhysicsList);

  // set mandatory user action class
  runManager->SetUserAction(new PrimaryGeneratorAction);

  // set MCTruth user action classes
  runManager->SetUserAction(new MCTruthTrackingAction);
  runManager->SetUserAction(new MCTruthEventAction);

  // Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/run/verbose 1");
  UI->ApplyCommand("/event/verbose 1");
  UI->ApplyCommand("/tracking/verbose 1");

  // configure MCTruth handling
  MCTruthConfig* config = new MCTruthConfig;
  config->SetMinE(1000.0);
  config->AddParticleType(11);
  MCTruthManager::GetInstance()->SetConfig(config);

  // start a run
  int numberOfEvent = 1;
  runManager->BeamOn(numberOfEvent);

  // job termination
  delete runManager;
  return 0;
}


