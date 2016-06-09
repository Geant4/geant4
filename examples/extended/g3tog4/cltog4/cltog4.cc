//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: cltog4.cc,v 1.6 2006/06/29 17:20:32 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// 
#include "G4ios.hh"
#include <fstream>
#include <cmath>

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

int main(int argc, char** argv)
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
    std::ifstream in(inFile);
    if (!in) {
      G4cout << "Cannot open input file """ << inFile << """" << G4endl;
      return 1;
    }
  }
  if (argc >= 3) {
    macroFile = argv[2];
    std::ifstream mac(macroFile);
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
