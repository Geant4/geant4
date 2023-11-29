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
#include <ctime>

#include "G4RunManagerFactory.hh"
#include "G4String.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"


int main(int argc, char **argv) {
  G4UIExecutive* pUi = nullptr;
  if (argc == 1) {
    G4cout << "WARNING! Please specify macro to be used: vis or run, or "
              "specify path to macro file. See README for more detail.\n";
    return 0;
  }
  if (argc > 1 && strcmp(argv[1], "vis") == 0) { pUi = new G4UIExecutive(argc, argv); }

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::MTwistEngine);
  // set random seed with system time
  G4Random::setTheSeed(time(nullptr));
  auto runManager= G4RunManagerFactory::CreateRunManager();
  runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());

  auto detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysicsList());
  runManager->SetUserInitialization(new ActionInitialization(detector));

  auto UImanager = G4UImanager::GetUIpointer();
  if (strcmp(argv[1], "run") == 0) {
    UImanager->ApplyCommand("/control/execute run.mac");
  }

  else if (strcmp(argv[1], "vis") == 0) {
    auto visManager = new G4VisExecutive;
    visManager->Initialize();
    UImanager->ApplyCommand("/control/execute vis.mac");
    pUi->SessionStart();
  }else {
    G4String filename = argv[1];
    UImanager->ApplyCommand("/control/execute " + filename);
  }
  delete runManager;
}
