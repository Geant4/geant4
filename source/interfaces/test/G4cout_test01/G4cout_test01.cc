//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//  G4cout_test02
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIsession.hh"

#include "ExN01DetectorConstruction.hh"
#include "ExN01PhysicsList.hh"
#include "ExN01PrimaryGeneratorAction.hh"
#include "MySession.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

int main()
{
  
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  runManager->SetUserInitialization(new ExN01DetectorConstruction);
  runManager->SetUserInitialization(new ExN01PhysicsList);

  // set mandatory user action class
  runManager->SetUserAction(new ExN01PrimaryGeneratorAction);

  // Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  // construct two sessions which receive G4cout/G4cerr


  G4UIsession * termSession = new G4UIterminal(new G4UItcsh); 


  // get the pointer to the UI manager and set verbosities
  UI->ApplyCommand("/run/verbose 1");
  UI->ApplyCommand("/event/verbose 1");
  UI->ApplyCommand("/tracking/verbose 1");

  // First, start tcsh session to see if all works well

  UI->SetCoutDestination(termSession);
  termSession->SessionStart();

  //  OK now end this interactive session
  delete termSession;
  //  delete termSession;

  G4cout << "term session end. now mySession starts to take a log file" << G4endl;

  // outputs are directed to the loggingSession from now on
  MySession * LoggingSession = new MySession;
  UI->SetCoutDestination(LoggingSession);
  // session->SessionStart(); // not required in this case
  // because we start a batch run directly here
  int numberOfEvent = 100;
  runManager->BeamOn(numberOfEvent);

  G4cout << "loggingSession is now deleted" << G4endl;

  delete LoggingSession;
  //  delete termSession;

  G4cout << "oooooooooooooooooooooooooooooooooooooooooooooooooo"<< G4endl;
  G4cout << "The output of G4cout/G4cerr is logged to a file" << G4endl;
    //       << LoggingSession->logFileName << "\"." << G4endl;
  G4cout << "Please have a look." << G4endl;
  G4cout << "oooooooooooooooooooooooooooooooooooooooooooooooooo" << G4endl;



  // job termination
  delete runManager;


  return 0;
}


