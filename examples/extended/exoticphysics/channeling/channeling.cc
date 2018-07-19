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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4ScoringManager.hh"
#include "G4UImanager.hh"

#ifdef G4MULTITHREADED
#include "UserActionInitialization.hh"
#else
#include "PrimaryGeneratorAction.hh"
#include "StackingAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#endif

#include "DetectorConstruction.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "FTFP_BERT.hh"

#include "G4EmStandardPhysics_option4_channeling.hh"
#include "G4ChannelingPhysics.hh"
#include "G4GenericBiasingPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
    // Detect interactive mode (if no arguments) and define UI session
    //
    G4UIExecutive* ui = 0;
    if ( argc == 1 ) {
      ui = new G4UIExecutive(argc, argv);
    }

    // Construct the default run manager
#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores() - 2);
#else
    G4RunManager* runManager = new G4RunManager;
#endif

    // Activate UI-command base scorer
    G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
    scManager->SetVerboseLevel(0);

    // Set mandatory initialization classes
    G4VModularPhysicsList* physlist= new FTFP_BERT();
    G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
    physlist->RegisterPhysics(new G4ChannelingPhysics());
    physlist->ReplacePhysics(new G4EmStandardPhysics_option4_channeling());
    //physlist->ReplacePhysics(new G4EmStandardPhysicsSS_channeling());
    biasingPhysics->PhysicsBiasAllCharged();
    physlist->RegisterPhysics(biasingPhysics);
    runManager->SetUserInitialization(physlist);

#ifndef G4MULTITHREADED
    // Set user action classes
    runManager->SetUserAction(new PrimaryGeneratorAction());
    runManager->SetUserAction(new EventAction());
    runManager->SetUserAction(new StackingAction());
    runManager->SetUserAction(new RunAction());
#else
    runManager->SetUserInitialization(new UserActionInitialization());
#endif

    runManager->SetUserInitialization(new DetectorConstruction());

    // Initialize visualization
    //
    G4VisManager* visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // Process macro or start UI session
    //
    if ( ! ui ) {
      // batch mode
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
    else {
      // interactive mode
      UImanager->ApplyCommand("/control/execute init_vis.mac");
      ui->SessionStart();
      delete ui;
    }


    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !

    delete visManager;
    delete runManager;

    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
