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

#include "eRositaActionInitialization.hh"
#include "eRositaDetectorConstruction.hh"
#include "eRositaEventAction.hh"
#include "eRositaPhysicsList.hh"
#include "eRositaPrimaryGeneratorAction.hh"
#include "eRositaRunAction.hh"
#include "eRositaSteppingAction.hh"
#include "eRositaSteppingVerbose.hh"

#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

auto main(int argumentCount, char** argumentVector) -> int
{
    // User verbose output class
    //
    auto *verbosity = new eRositaSteppingVerbose();
    G4VSteppingVerbose::SetInstance(verbosity);

    // Run manager
    //
    auto *runManager = G4RunManagerFactory::CreateRunManager();

    constexpr auto NUMBER_OF_THREADS{4};
    runManager->SetNumberOfThreads(NUMBER_OF_THREADS);

    G4cout << "***********************" << G4endl;
    G4cout << "** Seed: " << G4Random::getTheSeed() << " **" << G4endl;
    G4cout << "***********************" << G4endl;

    runManager->SetUserInitialization(new eRositaDetectorConstruction);
    runManager->SetUserInitialization(new eRositaPhysicsList);
    runManager->SetUserInitialization(new eRositaActionInitialization());

    // Initialize GEANT4 kernel
    //
    runManager->Initialize();

    // Get the pointer to the User Interface (UI) manager
    //
    auto *userInterfaceManager = G4UImanager::GetUIpointer();

    if (argumentCount != 1) { // batch mode
        G4String command = "/control/execute ";
        G4String fileName = argumentVector[1];
        userInterfaceManager->ApplyCommand(command + fileName);
    } else { // interactive mode : define visualization and UI terminal
        auto *visualizationManager = new G4VisExecutive;
        visualizationManager->Initialize();

        auto *userInterface = new G4UIExecutive(argumentCount, argumentVector);
        userInterfaceManager->ApplyCommand("/control/execute vis.mac");
        userInterface->SessionStart();
        
        delete userInterface;
        delete visualizationManager;
    }

    // Free the store:
    // user actions, physics list and detector description are
    // owned and deleted by the run manager,
    // so they should not be deleted in the main() program!

    delete runManager;
    delete verbosity;

    return 0;
}
