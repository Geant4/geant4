// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: cltog4.cc,v 1.1 2000-07-24 11:23:46 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "G4ios.hh"
#include "g4std/fstream"
#include <math.h>

// example header file includes

#include "G3toG4DetectorConstruction.hh"
#include "G3toG4RunAction.hh"
#include "G3toG4PrimaryGeneratorAction.hh"
#include "G3toG4PhysicsList.hh"
#include "G4LogicalVolume.hh"
#include "G3VolTable.hh"

// geant4 includes

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

G4int main(int argc, char** argv)
{
  G4String inFile;
  G4String macroFile = "";
  
  if (argc < 2) {
    G4cout << "Correct syntax:" << argv[0] << " <call_list_file> [ <macro_file> ]"
	   << G4endl;
    G4cout << "If only one argument is specified, macro file " << macroFile
	   << "will be used." << G4endl 
	   << "The second argument is used to override the default macro file"
	   << " name." << G4endl;
    
    return 1;
  }
  if (argc >= 2) {
    // Process the command line
    inFile = argv[1];
    G4std::ifstream in(inFile);
    if (!in) {
      G4cout << "Cannot open input file """ << inFile << """" << G4endl;
      return 1;
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
    G4cout << "Too many command line arguments (" << argc <<")" << G4endl;
    return 1;
  }
  // run manager
  G4RunManager * RunManager = new G4RunManager;
  // user initialization classes
  RunManager->SetUserInitialization(new G3toG4DetectorConstruction(inFile));
  RunManager->SetUserInitialization(new G3toG4PhysicsList);
  
  // user action classes
  RunManager->SetUserAction(new G3toG4RunAction);
  RunManager->SetUserAction(new G3toG4PrimaryGeneratorAction);

  G4UImanager * UI = G4UImanager::GetUIpointer();
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
  delete RunManager;
  return EXIT_SUCCESS;
}
