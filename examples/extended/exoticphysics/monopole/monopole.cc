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
/// \file exoticphysics/monopole/monopole.cc
/// \brief Main program of the exoticphysics/monopole example
//
// $Id: monopole.cc 110255 2018-05-17 14:15:38Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "DetectorConstruction.hh"
#include "G4MonopolePhysics.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "ActionInitialization.hh"

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::MixMaxRng);

  // Construct the default run manager
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(2);
  G4cout << "===== Monopole is started with "
         <<  runManager->GetNumberOfThreads() << " threads =====" << G4endl;
#else
  G4RunManager* runManager = new G4RunManager();
#endif

  //create physicsList
  // Physics List is defined via environment variable PHYSLIST
  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = factory.ReferencePhysList();

  // monopole physics is added
  G4MonopolePhysics * theMonopole = new G4MonopolePhysics();
  phys->RegisterPhysics(theMonopole);
  runManager->SetUserInitialization(phys);

  // visualization manager
  G4VisManager* visManager = new G4VisExecutive();
  visManager->Initialize();

  //get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Setup monopole
  G4String s = "";
  if(argc > 2) { s = argv[2]; }
  UImanager->ApplyCommand("/control/verbose 1");
  UImanager->ApplyCommand("/monopole/setup "+s);

  // set detector construction
  DetectorConstruction* det = new DetectorConstruction();
  runManager->SetUserInitialization(det);

  // set user action classes
  runManager->SetUserInitialization(new ActionInitialization(det));

  runManager->Initialize();

  if (ui)   // Define UI terminal for interactive mode
    {
      UImanager->ApplyCommand("/control/execute vis.mac");
      ui->SessionStart();
      delete ui;
    }
  else if (argc>1) // Batch mode with 1 or more files
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

  delete visManager;

  // job termination
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
