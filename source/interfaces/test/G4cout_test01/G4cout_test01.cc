//  G4cout_test01
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIsession.hh"

#include "ExN01DetectorConstruction.hh"
#include "ExN01PhysicsList.hh"
#include "ExN01PrimaryGeneratorAction.hh"
#include "MySession.hh"

#include "g4templates.hh"

int main()
{
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  // construct a session which receives G4cout/G4cerr

  MySession * LoggingSession = new MySession;
  UI->SetCoutDestination(LoggingSession);

  // session->SessionStart(); // not required in this case
  
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  runManager->SetUserInitialization(new ExN01DetectorConstruction);
  runManager->SetUserInitialization(new ExN01PhysicsList);

  // set mandatory user action class
  runManager->SetUserAction(new ExN01PrimaryGeneratorAction);

  // Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the UI manager and set verbosities
  UI->ApplyCommand("/run/verbose 1");
  UI->ApplyCommand("/event/verbose 1");
  UI->ApplyCommand("/tracking/verbose 1");

  // start a run
  int numberOfEvent = 3;
  runManager->BeamOn(numberOfEvent);


  // job termination
  delete runManager;
  cout << "oooooooooooooooooooooooooooooooooooooooooooooooooo"<< endl;
  cout << "The output of G4cout/G4cerr is logged to \""
       << LoggingSession->logFileName << "\"." << endl;
  cout << "Please have a look." << endl;
  cout << "oooooooooooooooooooooooooooooooooooooooooooooooooo" << endl;

  delete LoggingSession;
  return 0;
}


