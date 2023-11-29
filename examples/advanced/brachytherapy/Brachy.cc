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
// --------------------------------------------------------------
//                 GEANT 4 - Brachytherapy example
// --------------------------------------------------------------
//
// Code developed currently by:
//  S.Guatelli & D. Cutajar

//
//    *******************************
//    *                             *
//    *    Brachy.cc                *
//    *                             *
//    *******************************
//
//
#include "G4Types.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "BrachyActionInitialization.hh"
#include "G4VisExecutive.hh"
#include "BrachyDetectorConstruction.hh"
#include "BrachyPhysicsList.hh"
#include "BrachyPrimaryGeneratorAction.hh"
#include "G4SDManager.hh"
#include "Randomize.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"

#include "G4ScoringManager.hh"
#include "G4UIExecutive.hh"

#include "G4ScoringManager.hh"
#include "BrachyUserScoreWriter.hh"

int main(int argc ,char ** argv)

{
 // Construct the default run manager
 //
  auto* pRunManager = G4RunManagerFactory::CreateRunManager();
  G4int nThreads = 4;
  pRunManager->SetNumberOfThreads(nThreads);
 
  G4cout << "***********************" << G4endl;
  G4cout << "*** Seed: " << G4Random::getTheSeed() << " ***" << G4endl;
  G4cout << "***********************" << G4endl;
 
  // Access to the Scoring Manager pointer

  G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManager();

  // Overwrite the default output file with user-defined one
  scoringManager->SetScoreWriter(new BrachyUserScoreWriter());

  // Initialize the physics component
  pRunManager -> SetUserInitialization(new BrachyPhysicsList);

  // Initialize the detector component
  auto pDetectorConstruction = new  BrachyDetectorConstruction();
  pRunManager -> SetUserInitialization(pDetectorConstruction);

  // User action initialization

 auto actions = new BrachyActionInitialization();
  pRunManager->SetUserInitialization(actions);

  //Initialize G4 kernel
 // pRunManager -> Initialize();

  // Visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if (argc == 1)   // Define UI session for interactive mode.
    {
      auto ui = new G4UIExecutive(argc, argv);
      G4cout << " UI session starts ..." << G4endl;
      UImanager -> ApplyCommand("/control/execute VisualisationMacro.mac");
      ui -> SessionStart();
      delete ui;
    }
  else           // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager -> ApplyCommand(command+fileName);
    }

  // Job termination
 // Close the root file

  delete visManager;
  delete pRunManager;

  return 0;
}
