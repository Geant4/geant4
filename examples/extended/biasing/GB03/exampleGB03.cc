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

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4UImanager.hh"

#include "GB03DetectorConstruction.hh"
#include "GB03PrimaryGeneratorAction.hh"
#include "GB03ActionInitialization.hh"
#include "FTFP_BERT.hh"
#include "G4GenericBiasingPhysics.hh"

int main(int argc,char** argv)
{
  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // -- Construct the run manager : MT or sequential one
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  G4cout << "      ********** Run Manager constructed in MT mode ************ " << G4endl;
  // -- Choose 4 threads:
  runManager->SetNumberOfThreads(4);
#else
  G4RunManager * runManager = new G4RunManager;
  G4cout << "      ********** Run Manager constructed in sequential mode ************ " << G4endl;
#endif

  // -- Set mandatory initialization classes
  G4VUserDetectorConstruction* detector = new GB03DetectorConstruction;
  runManager->SetUserInitialization(detector);
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  // -- Setup biasing for neutron only. As will not bias any
  // -- physics process, require the G4GenericBiasingPhysics()
  // -- to modify the physics list just for "non-physics" biasing:
  G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
  biasingPhysics->NonPhysicsBias("neutron");
  physicsList->RegisterPhysics(biasingPhysics);
  runManager->SetUserInitialization(physicsList);

  // -- Set user action initialization class
  G4VUserActionInitialization* actions = new GB03ActionInitialization;
  runManager->SetUserInitialization(actions);

  // Visualization manager
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // Initialize G4 kernel
  //
  runManager->Initialize();

  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (ui)   // Define UI session for interactive mode
    {
      UImanager->ApplyCommand("/control/execute vis.mac");
      ui->SessionStart();
      delete ui;
    }
  else           // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

  delete visManager;
  delete runManager;

  return 0;
}
