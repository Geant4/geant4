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
/// \file field/field04/field04.cc
/// \brief Main program of the field/field04 example
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef WIN32
#  include <unistd.h>
#endif

#include "G4Types.hh"

#include "G4RunManagerFactory.hh"
#include "F04SteppingVerbose.hh"

#include "F04ActionInitialization.hh"
#include "F04DetectorConstruction.hh"
#include "F04PhysicsList.hh"

#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace
{
void PrintUsage()
{
  G4cerr << " Usage: " << G4endl;
  G4cerr << " field04 [-m macro ] [-p physicsList] [-r randomSeed] [-s preinit|idle]" << G4endl;
}
}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // Evaluate arguments
  // argc holds the number of arguments (including the name) on the command line
  // -> it is ONE when only the name is  given !!!
  // argv[0] is always the name of the program
  // argv[1] points to the first argument, and so on
  //
  if (argc > 7) {
    PrintUsage();
    return 1;
  }

  G4String macro;
  G4String physicsList = "QGSP_BERT";
  G4int randomSeed = 1234;
  G4String startPhase = "idle";
  for (G4int i = 1; i < argc; i = i + 2) {
    if (G4String(argv[i]) == "-m")
      macro = argv[i + 1];
    else if (G4String(argv[i]) == "-r")
      randomSeed = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-p")
      physicsList = argv[i + 1];
    else if (G4String(argv[i]) == "-s")
      startPhase = argv[i + 1];
    else {
      PrintUsage();
      return 1;
    }
  }

  // Instantiate G4UIExecutive if there are no arguments (interactive mode)
  //
  G4UIExecutive* ui = nullptr;
  if (!macro.size()) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Setting the application-specific SteppingVerbose
  //
  auto verbosity = new F04SteppingVerbose;

  // Construct the default run manager
  //
  auto runManager = G4RunManagerFactory::CreateRunManager();

  G4Random::setTheSeed(randomSeed);

  // Set mandatory initialization classes
  //
  // Detector construction
  auto detector = new F04DetectorConstruction();
  runManager->SetUserInitialization(detector);
  // Physics list
  runManager->SetUserInitialization(new F04PhysicsList(physicsList));
  // User action initialization
  runManager->SetUserInitialization(new F04ActionInitialization(detector));

  // Initialize G4 kernel
  //
  // runManager->Initialize();

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if (macro.size()) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + macro);
  }
  else {
    // interactive mode : define UI session
    if (startPhase == "preinit") {
      // start in PreInit> phase if requested
      G4cout << "At the prompt, issue commands to set up detector & field, then:" << G4endl;
      G4cout << "/run/initialize" << G4endl;
      G4cout << "Then if you want a viewer:" << G4endl;
      G4cout << "/control/execute vis.mac" << G4endl;
      G4cout << "Then: " << G4endl;
      G4cout << "/run/beamOn … etc." << G4endl;
    }
    else {
      // perform initialization and draw geometry
      UImanager->ApplyCommand("/control/execute init_vis.mac");
    }
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  delete verbosity;
  delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
