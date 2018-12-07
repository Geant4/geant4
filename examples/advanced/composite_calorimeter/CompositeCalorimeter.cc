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
///////////////////////////////////////////////////////////////////////////////
// File: CompositeCalorimeter.cc
// Description: Main function for Geant4 application HCAL Test-BEAM H2-96
///////////////////////////////////////////////////////////////////////////////

#include "CCalDetectorConstruction.hh"
#include "CCalActionInitializer.hh"
#include "CCalAnalysis.hh"

#include "G4PhysListFactory.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4VisExecutive.hh"

#include "G4UIExecutive.hh"

int main(int argc,char** argv) {

  G4VisManager *visManager = nullptr;

  G4RunManager * runManager = new G4RunManager;
  runManager->SetUserInitialization(new CCalDetectorConstruction);

  G4PhysListFactory factory;
  runManager->SetUserInitialization(factory.ReferencePhysList());

  runManager->SetUserInitialization(new CCalActionInitializer());
  
  G4UImanager * UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/CCal/generator/verbose 2");
  UImanager->ApplyCommand("/gun/position -1380. 0. 0. mm");
  UImanager->ApplyCommand("/gun/direction 1. 0. 0.");
  UImanager->ApplyCommand("/gun/energy 100 GeV");

  // Define (G)UI terminal for interactive mode
  if (argc==1) {  // No arguments - interactive assumed.

    visManager = new G4VisExecutive;
    visManager->Initialize();     

    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    if (ui->IsGUI()) {
      // Customize the menubar with a macro file :
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    G4cout <<" Run initializing ..."<<G4endl;
    UImanager->ApplyCommand("/process/verbose 0");
    UImanager->ApplyCommand("/run/verbose 2");
    UImanager->ApplyCommand("/run/initialize");

    // Create empty scene
    G4String visCommand = "/vis/scene/create";
    UImanager->ApplyCommand(visCommand);
    
    // Choose one default viewer 
    // (the user can always change it later on) 
    // visCommand = "/vis/open DAWNFILE";
    // visCommand = "/vis/open VRML2FILE";
    visCommand = "/vis/open OGL";
    UImanager->ApplyCommand(visCommand);
    
    visCommand = "/vis/viewer/flush";
    UImanager->ApplyCommand(visCommand);
    visCommand = "/tracking/storeTrajectory 1";
    UImanager->ApplyCommand(visCommand);

    G4cout <<"Now, please, apply beamOn command..."<<G4endl;
    ui->SessionStart();
    delete ui;

  } else {

    // Batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  //Close-out analysis:
  // Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
  // Complete clean-up
  delete G4AnalysisManager::Instance();
  delete runManager;
  delete visManager;
  return 0;
}
