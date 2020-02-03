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
/// \file persistency/gdml/G03/gdml_ext.cc
/// \brief Main program of the persistency/gdml/G03 example
//
//
<<<<<<< HEAD
// $Id: gdml_ext.cc 68025 2013-03-13 13:43:46Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
//
// --------------------------------------------------------------
//      GEANT 4 - read_ext
//
// --------------------------------------------------------------

// Geant4 includes
//
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "globals.hh"

// A pre-built physics list
//
#include "QGSP_BERT.hh"

// Example includes
//
#include "G03DetectorConstruction.hh"
#include "G03PrimaryGeneratorAction.hh"
#include "G03RunAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc, char** argv)
{
       
  // Construct the default run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // Set mandatory initialization and user action classes
  //
  G03DetectorConstruction* detector = new G03DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new QGSP_BERT);
  runManager->SetUserAction(new G03PrimaryGeneratorAction);
  G03RunAction* runAction = new G03RunAction;
  runManager->SetUserAction(runAction);
      
  // Initialisation of runManager via macro for the interactive mode
  // This gives possibility to give different names for GDML file to READ
 
#ifdef G4VIS_USE
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // run initialisation macro
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if ( argc==1 )   // Define UI session for interactive mode. 
  {
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute vis.mac");     
#endif
    ui->SessionStart();
    delete ui;
#endif
  }
  else             // Batch mode
  { 
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    G4String command = "/control/execute "; 
    G4String fileName = argv[1]; 
    UImanager->ApplyCommand(command+fileName); 
    ui->SessionStart();
    delete ui;
#endif
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  
  // Job termination
  //
  delete runManager;

  return 0;
}
