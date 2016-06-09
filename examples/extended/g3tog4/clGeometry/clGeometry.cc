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
/// \file g3tog4/clGeometry/clGeometry.cc
/// \brief Main program of the g3tog4/clGeometry example
//
//
// $Id$
//
// 

// package includes
#include "G3toG4DetectorConstruction.hh"

// common package includes
#include "ExG4EventAction01.hh"
#include "ExG4RunAction01.hh"
#include "ExG4PrimaryGeneratorAction01.hh"

// geant4 includes
#include "FTFP_BERT.hh"
#include "G3VolTable.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

// visualization
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

// (G)UI
#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  G4String inFile;
  G4String macroFile = "";
  G4bool batchMode = false;
    
  if (argc < 2 || argc >= 4) {
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
    G4cout << "Geometry data file: " << inFile << G4endl;
    std::ifstream in(inFile);
    if (!in) {
      G4cerr << "Cannot open input file \"" << inFile << "\"" << G4endl;
      return EXIT_FAILURE;
    }
  }
  if (argc == 3) {
    batchMode = true;
    macroFile = argv[2];
    std::ifstream mac(macroFile);
    if (!mac) {
      G4cout << "Cannot open macro file """ << macroFile << """" << G4endl;
      return 2;
    }
  }
    
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;
    
  // set mandatory initialization classes
  runManager->SetUserInitialization(new G3toG4DetectorConstruction(inFile));
  runManager->SetUserInitialization(new FTFP_BERT);
    
  // set user action classes
  runManager->SetUserAction(new ExG4PrimaryGeneratorAction01());
  runManager->SetUserAction(new ExG4RunAction01);
  runManager->SetUserAction(new ExG4EventAction01);
    
  // Initialize G4 kernel
  //
  runManager->Initialize();

  //----------------
  // Visualization:
  //----------------

  // the pointer to the User Interface manager 
  G4UImanager* uiManager = G4UImanager::GetUIpointer();  

  if ( batchMode ) {
    // batch mode
    uiManager->ApplyCommand("/control/execute init.mac"); 
    G4String command = "/control/execute ";
    uiManager->ApplyCommand(command+macroFile);
  }
  else {
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager -> Initialize();
    uiManager->ApplyCommand("/control/execute init_vis.mac");
#else
    uiManager->ApplyCommand("/control/execute init.mac"); 
#endif
    ui->SessionStart();
    delete ui;
#endif

#ifdef G4VIS_USE
    delete visManager;
#endif
  } 

  delete runManager;
  return EXIT_SUCCESS;
}
