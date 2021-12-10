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
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo
//
//   *******************************************************
//   *                  Ultra.cc
//   *******************************************************

#include "G4Types.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "UltraActionInitializer.hh"
#include "UltraDetectorConstruction.hh"
#include "UltraPhysicsList.hh"


int main(int argc,char** argv) {

//choose the Random engine from CLHEP
//(lets use C++ implementation of Jame's RANLUX generator)

  G4Random::setTheEngine(new CLHEP::RanluxEngine);

  auto* runManager = G4RunManagerFactory::CreateRunManager();
  G4int nThreads = 4;
  runManager->SetNumberOfThreads(nThreads);

  // UserInitialization classes - mandatory
  runManager->SetUserInitialization(new UltraDetectorConstruction);
  runManager->SetUserInitialization(new UltraPhysicsList());

  // UserAction classes - optional
  runManager->SetUserInitialization(new UltraActionInitializer());

  // Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Initialise visualization
  auto* visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the Pointer to the UI Manager
  auto* UImanager = G4UImanager::GetUIpointer();

  // User interactions
  // Define (G)UI for interactive mode
  if(argc==1)
  {
    UImanager->ApplyCommand("/control/execute VisDirect.mac");
    ui->SessionStart();
    delete ui;
  }
  else  // Batch mode
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  delete visManager;
  delete runManager;

  return 0;
}
