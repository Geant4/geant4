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
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc, argv);

  // Use SteppingVerbose with Unit
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  if (argc == 2) {
    G4Exception("argument", "Fatal error in Argument", FatalErrorInArgument,
                "Type of simulation must be precised: either ./stim_pixe_tomography -p "
                "pixe3d.mac or ./stim_pixe_tomography -s pixe3d.mac");
  }
  G4bool isPIXE = false;
  if (argc > 2) {
    G4String s;
    s = argv[1];
    if (s == "-p")
      isPIXE = true;
    else if (s == "-s")
      isPIXE = false;
    else {
      G4Exception("argument", "Fatal error in Argument", FatalErrorInArgument,
                  "Type of simulation must be precised: /stim_pixe_tomography -p "
                  "pixe3d.mac or ./stim_pixe_tomography -s pixe3d.mac");
    }
  }
  G4int nThreads = 4;
  if (argc == 4) nThreads = G4UIcommand::ConvertToInt(argv[3]);
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  runManager->SetNumberOfThreads(nThreads);

  // UserInitialization classes - mandatory
  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new PhysicsList());
  runManager->SetUserInitialization(new ActionInitialization(isPIXE));

  // initialize visualization
  auto visManager = new G4VisExecutive;
  visManager->Initialize();

  // get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // User interactions
  // Define (G)UI for interactive mode
  if (ui) {
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }
  else {
    G4String command;
    command = "/control/execute ";
    G4String fileName;
    fileName = argv[2];
    UImanager->ApplyCommand(command + fileName);
  }

  // job termination
  delete visManager;
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
