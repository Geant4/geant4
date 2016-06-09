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
/// \file visualization/userVisAction/userVisAction.cc
/// \brief Main program of the visualization/userVisAction example
//
//
// $Id$
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B1DetectorConstruction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1RunAction.hh"
#include "B1EventAction.hh"
#include "B1SteppingAction.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "QGSP_BIC_EMY.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#include "UVA_VisAction.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes (mandatory)
  runManager->SetUserInitialization(new B1DetectorConstruction);
  runManager->SetUserInitialization(new QGSP_BIC_EMY);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  // Set User Vis Action...
  visManager->SetUserAction
    (new UVA_VisAction, G4VisExtent(-20*cm,-10*cm,-25*cm,-15*cm,20*cm,40*cm));
#endif
   
  // UserAction classes
  runManager->SetUserAction(new B1PrimaryGeneratorAction);
  runManager->SetUserAction(new B1RunAction);  
  runManager->SetUserAction(new B1EventAction);
  runManager->SetUserAction(new B1SteppingAction);

  //Initialize G4 kernel
  runManager->Initialize();
      
  //get the pointer to the User Interface manager 
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  

  if(argc==1)
  // Define (G)UI terminal for interactive mode  
  { 
#ifdef G4UI_USE
    G4UIExecutive * ui = new G4UIExecutive(argc, argv);
    UImanager->ApplyCommand("/control/execute vis.mac");    
    ui->SessionStart();
    delete ui;
#endif
  }
  else
  // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

