// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: clGeometry.cc,v 1.1 2000-07-24 11:23:40 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

// controls whether drawing is to be Done or not

#include "g4std/fstream"
#include <math.h>
#include "G4ios.hh"

// package includes

#include "G3toG4DetectorConstruction.hh"
#include "G3toG4RunAction.hh"
#include "G3toG4PrimaryGeneratorAction.hh"
#include "G3toG4PhysicsList.hh"
#include "G3toG4EventAction.hh"
#include "G4LogicalVolume.hh"
#include "G3VolTable.hh"

// geant4 includes

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

// visualization
#ifdef G4VIS_USE
#include "G3toG4VisManager.hh"
#endif

G4int main(int argc, char** argv)
{
  G4String inFile;
  G4String macroFile = "";
    
  if (argc < 2) {
    G4cerr << "clGeometry: Correct syntax: clGeometry <call_list_file> [ <macro_file> ]"
	   << G4endl;
    G4cerr << "If only one argument is specified, interactive mode will be "
	   << "entered." << G4endl << "The second argument, if specified, is "
	   << "the name of the macro file (batch mode)." << G4endl;
        
    return EXIT_FAILURE;
  }
  if (argc >= 2) {
    // Process the command line
    inFile = argv[1];
    G4std::ifstream in(inFile);
    if (!in) {
      G4cerr << "Cannot open input file \"" << inFile << "\"" << G4endl;
      return EXIT_FAILURE;
    }
  }
  if (argc >= 3) {
    macroFile = argv[2];
    G4std::ifstream mac(macroFile);
    if (!mac) {
      G4cout << "Cannot open macro file """ << macroFile << """" << G4endl;
      return 2;
    }
  }
  if (argc >= 4) {
    G4cerr << "Too many command line arguments (" << argc <<")" << G4endl;
    return EXIT_FAILURE;
  }
    
  // Construct the default run manager
  G4RunManager* RunManager = new G4RunManager;
    
  // set mandatory initialization classes
  RunManager->SetUserInitialization(new G3toG4DetectorConstruction(inFile));

  G3toG4PhysicsList* thePhysicsList = new G3toG4PhysicsList;

  // set verbosity of PhysicsList
  thePhysicsList->SetVerboseLevel(2);
  RunManager->SetUserInitialization(thePhysicsList);
    
  //----------------
  // Visualization:
  //----------------

#ifdef G4VIS_USE
  G4VisManager* VisManager = new G3toG4VisManager;
  VisManager -> Initialize();
#endif
    
  // set user action classes

  RunManager->SetUserAction(new G3toG4RunAction);

  G3toG4EventAction* theEventAction = new G3toG4EventAction;
  theEventAction->SetDrawFlag("all");
  RunManager->SetUserAction(theEventAction);

  RunManager->SetUserAction(new G3toG4PrimaryGeneratorAction);
    
  // the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  // set some additional defaults and initial actions
    
  UI->ApplyCommand("/control/verbose 1");
  UI->ApplyCommand("/run/verbose 1");
  UI->ApplyCommand("/tracking/verbose 1");
  UI->ApplyCommand("/tracking/storeTrajectory 1");
  UI->ApplyCommand("/run/initialize");

  G4bool batch_mode = macroFile != "";
    
  if(!batch_mode) {
    G4UIsession * session = new G4UIterminal;
    if (session != 0) {
      session->SessionStart();
      delete session;
      //      G4cout << "deleted G4UITerminal..." << G4endl;
    }
  }
  else {
    // Batch mode
    G4String command = "/control/execute ";
    UI->ApplyCommand(command+macroFile);
  }
#ifdef G4VIS_USE
  if (VisManager !=0) delete VisManager;
#endif
  delete RunManager;
  return EXIT_SUCCESS;
}

























