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
/// \file microtrack.cc
/// \brief Main program of the dna/microtrack example

// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 45, (2018) e722-e739
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// Contact author: M.A. Cortes-Giraldo, Universidad de Sevilla, Spain
//                 (miancortes@us.es)
// Other authors: M. Galocha-Oliva, A. Baratto-Roldan
//
// Example ref. paper: A. Baratto-Roldan et al., Front Phys 9 (2021) 726787
//

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc, argv);

  // choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Use SteppingVerbose with Unit
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  // Creating run manager
  auto runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  if (argc == 3) {
    G4int nThreads = G4UIcommand::ConvertToInt(argv[2]);
    runManager->SetNumberOfThreads(nThreads);
  }

  // set mandatory initialization classes
  DetectorConstruction* det = new DetectorConstruction;
  runManager->SetUserInitialization(det);

  PhysicsList* phys = new PhysicsList;
  runManager->SetUserInitialization(phys);

  runManager->SetUserInitialization(new ActionInitialization(/*det, phys*/));

  // initialize visualization
  G4VisManager* visManager = nullptr;

  // get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (ui) {
    // interactive mode
    visManager = new G4VisExecutive;
    visManager->Initialize();
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }
  else {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }

  // job termination
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
