#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef G4VIS_USE
#include "GammaRayTelVisManager.hh"
#endif

#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelPhysicsList.hh"
#include "GammaRayTelPrimaryGeneratorAction.hh"
#include "GammaRayTelRunAction.hh"
#include "GammaRayTelEventAction.hh"

// This global file is used to store relevant data for
// analysis with a separate program
ofstream outFile("prova.dat");

int main()
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  GammaRayTelDetectorConstruction* detector = new GammaRayTelDetectorConstruction;
  runManager->SetUserInitialization(detector);

  runManager->SetUserInitialization(new GammaRayTelPhysicsList);

  G4UIsession* session=0;
  session = new G4UIterminal; // interactive

#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new GammaRayTelVisManager;
  visManager->Initialize();
#endif

  // set mandatory user action classes
  runManager->SetUserAction(new GammaRayTelPrimaryGeneratorAction(detector));
  runManager->SetUserAction(new GammaRayTelRunAction);

  // set optional user action classes
  GammaRayTelEventAction* eventAction = new GammaRayTelEventAction;
  runManager->SetUserAction(eventAction);

  // Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (session)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
      // prerunGammaRayTel.mac is loaded by default
      UI->ApplyCommand("/control/execute prerunGammaRayTel.mac");
    
      session->SessionStart();
      delete session;
    }

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  return 0;
}






