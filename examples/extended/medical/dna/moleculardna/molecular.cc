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
/// \file molecular.cc
/// \brief Molecular level simulation of DNA

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PhysicsList.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

#include <ctime>
#include "G4MoleculeGun.hh"
#include "G4DNAChemistryManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
    void PrintUsage() {
      G4cout << " Usage: " << G4endl;
      G4cout << " molecular [-m macro ] [-t nThreads] [-p PhysicsList]"
             << G4endl;
      G4cout << "   -p is the G4DNA Physics List option. Default (0) is"
             << " G4EmDNAPhysics" << G4endl;
      G4cout << "   note: -t option is available only for multi-threaded mode."
             << G4endl;
    }
}

int main(int argc, char **argv) {
  if (argc > 7) {
    PrintUsage();
    return 1;
  }

  G4String macro;
  G4int phys_option = 2;

  G4int nThreads = 2;
  for (G4int ii = 1; ii < argc; ii = ii + 2) {
    if (G4String(argv[ii]) == "-m") {
      macro = argv[ii + 1];
    } else if (G4String(argv[ii]) == "-p") {
      phys_option = G4UIcommand::ConvertToInt(argv[ii + 1]);
    } else if (G4String(argv[ii]) == "-t") {
      nThreads = G4UIcommand::ConvertToInt(argv[ii + 1]);
    } else {
      PrintUsage();
      return 1;
    }
  }

  G4UIExecutive *ui = nullptr;
  if (!macro.size()) {
    ui = new G4UIExecutive(argc, argv);
  }

  /*
  time_t timeStart;
  time(&timeStart);
  long seed = timeStart;
  G4Random::setTheSeed(seed);
  */

  auto* runManager = G4RunManagerFactory::CreateRunManager();
  if (nThreads > 0) {
    runManager->SetNumberOfThreads(nThreads);
  }

  runManager->SetUserInitialization(new DetectorConstruction());
  G4VModularPhysicsList *physicsList = new PhysicsList(phys_option);
  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new ActionInitialization());
  G4DNAChemistryManager::Instance()->Initialize();

  G4VisExecutive* visManager = nullptr;

  // Get the pointer to the User Interface manager
  G4UImanager *UImanager = G4UImanager::GetUIpointer();
  if (!macro.empty()) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + macro);
  } else {
    visManager = new G4VisExecutive;
    visManager->Initialize();
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
    delete visManager;
  }
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
