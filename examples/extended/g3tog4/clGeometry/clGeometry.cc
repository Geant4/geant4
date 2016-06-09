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
// $Id: clGeometry.cc,v 1.9 2010-11-04 21:47:43 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

// controls whether drawing is to be Done or not

#include <fstream>
#include <cmath>
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

// visualization
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

// (G)UI
#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc, char** argv)
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
    std::ifstream in(inFile);
    if (!in) {
      G4cerr << "Cannot open input file \"" << inFile << "\"" << G4endl;
      return EXIT_FAILURE;
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
  G4VisManager* VisManager = new G4VisExecutive;
  VisManager -> Initialize();
#endif
    
  // set user action classes

  RunManager->SetUserAction(new G3toG4RunAction);

  G3toG4EventAction* theEventAction = new G3toG4EventAction;
  theEventAction->SetDrawFlag("all");
  RunManager->SetUserAction(theEventAction);

  RunManager->SetUserAction(new G3toG4PrimaryGeneratorAction);
    
  // the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  // set some additional defaults and initial actions
    
  UImanager->ApplyCommand("/control/verbose 1");
  UImanager->ApplyCommand("/run/verbose 1");
  UImanager->ApplyCommand("/tracking/verbose 1");
  UImanager->ApplyCommand("/tracking/storeTrajectory 1");
  UImanager->ApplyCommand("/run/initialize");

  G4bool batch_mode = macroFile != "";
    
  if(!batch_mode) {
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute vis.mac");
#endif
    ui->SessionStart();
    delete ui;
#endif
  } else {
    // Batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macroFile);
  }
#ifdef G4VIS_USE
  delete VisManager;
#endif
  delete RunManager;
  return EXIT_SUCCESS;
}
