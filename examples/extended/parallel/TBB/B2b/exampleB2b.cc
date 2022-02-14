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
/// \file parallel/TBB/B2b/exampleB2b.cc
/// \brief Main program of the B2b example

#include "B2bDetectorConstruction.hh"
#include "B2ActionInitialization.hh"

#include "G4Threading.hh"
#include "G4TaskRunManager.hh"

#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // Instantiate G4UIExecutive if there are no arguments (interactive mode)
  G4UIExecutive* ui = 0;
  if(argc == 1)
  {
    ui = new G4UIExecutive(argc, argv);
  }

  // Choose the Random engine

  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  auto runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::TBB);

  //==== Geant4 specific stuff, from now up to END-G4 comment is copy from MT
  // example
  // Set mandatory initialization classes

  runManager->SetUserInitialization(new B2bDetectorConstruction());

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  // Set user action classes

  runManager->SetUserInitialization(new B2ActionInitialization());

  // Initialize G4 kernel

  runManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if(!ui)  // batch mode
  {
    G4String command  = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }
  else
  {
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init.mac");
    if(ui->IsGUI())
      UImanager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();
    delete ui;
  }

  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
