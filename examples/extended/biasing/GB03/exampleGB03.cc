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
/// \file exampleGB03.cc
/// \brief Main program of the biasing/GB03 example

#include "FTFP_BERT.hh"
#include "GB03ActionInitialization.hh"
#include "GB03DetectorConstruction.hh"
#include "GB03PrimaryGeneratorAction.hh"

#include "G4GenericBiasingPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace
{
void PrintUsage()
{
  G4cerr << " Usage: " << G4endl;
  G4cerr << " ./exampleGB04 [-m macro ] "
         << " [-b biasing {'on','off'}]"
         << "\n or\n ./exampleGB04 [macro.mac]" << G4endl;
}
}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // Evaluate arguments
  //
  if (argc > 5) {
    PrintUsage();
    return 1;
  }

  G4String macro("");
  G4String onOffBiasing("");
  if (argc == 2)
    macro = argv[1];
  else {
    for (G4int i = 1; i < argc; i = i + 2) {
      if (G4String(argv[i]) == "-m")
        macro = argv[i + 1];
      else if (G4String(argv[i]) == "-b")
        onOffBiasing = argv[i + 1];
      else {
        PrintUsage();
        return 1;
      }
    }
  }

  if (onOffBiasing == "") onOffBiasing = "on";

  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  if (macro == "") {
    ui = new G4UIExecutive(argc, argv);
  }

  // -- Construct the run manager : MT or sequential one
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  runManager->SetNumberOfThreads(4);

  G4bool biasingFlag = (onOffBiasing == "on");

  // -- Set mandatory initialization classes
  G4VUserDetectorConstruction* detector = new GB03DetectorConstruction(biasingFlag);
  runManager->SetUserInitialization(detector);
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  // -- Setup biasing for neutron only. As will not bias any
  // -- physics process, require the G4GenericBiasingPhysics()
  // -- to modify the physics list just for "non-physics" biasing:
  if (biasingFlag) {
    auto biasingPhysics = new G4GenericBiasingPhysics();
    biasingPhysics->NonPhysicsBias("neutron");
    physicsList->RegisterPhysics(biasingPhysics);
    G4cout << "      ********************************************************* " << G4endl;
    G4cout << "      ********** processes are wrapped for biasing ************ " << G4endl;
    G4cout << "      ********************************************************* " << G4endl;
  }
  else {
    G4cout << "      ************************************************* " << G4endl;
    G4cout << "      ********** processes are not wrapped ************ " << G4endl;
    G4cout << "      ************************************************* " << G4endl;
  }
  runManager->SetUserInitialization(physicsList);

  // -- Set user action initialization class
  G4VUserActionInitialization* actions = new GB03ActionInitialization;
  runManager->SetUserInitialization(actions);

  // Visualization manager
  //
  auto visManager = new G4VisExecutive;
  visManager->Initialize();

  // Initialize G4 kernel
  //
  runManager->Initialize();

  // Get the pointer to the User Interface manager
  //
  auto UImanager = G4UImanager::GetUIpointer();

  if (ui)  // Define UI session for interactive mode
  {
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  }
  else  // Batch mode
  {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + macro);
  }

  delete visManager;
  delete runManager;

  return 0;
}
