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
/// \file biasing/ReverseMC01/exampleRMC01.cc
/// \brief Main program of the biasing/ReverseMC01 example
//
//
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleRMC01
//
// --------------------------------------------------------------
// Comments
//
// This example intends to show how to use the adjoint/reverse simulation mode.
// 
// --------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4AdjointSimManager.hh"
#include "RMC01DetectorConstruction.hh"
#include "RMC01PrimaryGeneratorAction.hh"
#include "RMC01EventAction.hh"
#include "RMC01RunAction.hh"
#include "G4AdjointPhysicsList.hh"

#include "G4AdjointSimManager.hh"
#include "RMC01AdjointEventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {
   
  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Construct the default run manager
  G4RunManager * theRunManager = new G4RunManager;
   
  
  RMC01DetectorConstruction* detector = new RMC01DetectorConstruction();
  
  //Physics and geometry are declared as in a normal G4 application
  //--------------------------------------------
  theRunManager->SetUserInitialization(detector);
  theRunManager->SetUserInitialization(new G4AdjointPhysicsList);

  theRunManager->SetUserAction(new RMC01PrimaryGeneratorAction);
  theRunManager->SetUserAction(new RMC01EventAction);
  RMC01RunAction* theRunAction = new RMC01RunAction;
  theRunManager->SetUserAction(theRunAction);
 
  //The adjoint simulation manager will control the Reverse MC mode 
  //---------------------------------------------------------------
  
  G4AdjointSimManager* theAdjointSimManager
                            = G4AdjointSimManager::GetInstance();
  
  //It is possible to define action that will be used during
  //                                        the adjoint tracking phase
  
  theAdjointSimManager->SetAdjointRunAction(theRunAction);
  //theAdjointSimManager->SetAdjointEventAction(new RMC01AdjointEventAction);
  theAdjointSimManager->SetAdjointEventAction(new RMC01EventAction);
  
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if (!ui)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);    
    }
  else
    {  // interactive mode : define UI session
      ui->SessionStart();
      delete ui;
    }

  // job termination
  delete visManager;
  delete theRunManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
