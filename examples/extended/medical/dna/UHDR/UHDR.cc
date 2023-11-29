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
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "PhysicsList.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

int main(int argc, char **argv) {
  G4int seed(0);
  if (argc > 2) {
    seed = std::stoi(argv[2]);
    long enterseed = seed;
    G4Random::setTheSeed(enterseed);
    G4Random::showEngineStatus();
    G4cout << "seed number : " << seed << G4endl;
  }
  auto runManager = G4RunManagerFactory::CreateRunManager();
  // Set mandatory initialization classes

  auto pDetector = new DetectorConstruction();
  runManager->SetUserInitialization(pDetector);
  runManager->SetUserInitialization(new PhysicsList(pDetector));
  runManager->SetUserInitialization(new ActionInitialization(pDetector));

  ///
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetNtupleDirectoryName("ntuple");
  // open output file
  //
  std::string name = "Dose_" + std::to_string(seed);
  G4bool fileOpen = analysisManager->OpenFile(name);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " << name << G4endl;
    return 0;
  }
  analysisManager->SetFirstNtupleId(1);
  ///////
  G4UImanager *UI = G4UImanager::GetUIpointer();
  if (argc > 1) // batch mode
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command + fileName);
  } else {
    UI->ApplyCommand("/control/execute beam.in");
  }
  analysisManager->Write();
  analysisManager->CloseFile();
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
