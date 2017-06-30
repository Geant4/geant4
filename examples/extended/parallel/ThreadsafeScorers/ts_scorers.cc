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
/// \file parallel/ThreadsafeScorers/ts_scorers.cc
/// \brief Main of the ThreadsafeScorers example
//
//
// $Id: ts_scorers.cc 93110 2015-11-05 08:37:42Z jmadsen $
//
//
/// ts_scorers example shows how to use global scorers. The benefit of using
///     global scorers in memory-savings for problems with very large amounts
///     of scoring volumes. Additionally, the global scorers are more precise
///     w.r.t. the serial solution because of the lack of compounding
///     round-off error from multiple threads
///
/// In this example, the global scorers are implemented as static member
///     variables in TSRun because TSRun is thread-local. The G4atomic
///     class is the core of the thread-safe scorers and can be uses
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifdef G4MULTITHREADED
    #include "G4MTRunManager.hh"
    #include "G4Threading.hh"
    typedef G4MTRunManager RunManager;
#else
    #include "G4RunManager.hh"
    typedef G4RunManager RunManager;
#endif

#include "Randomize.hh"

// User Defined Classes
#include "TSDetectorConstruction.hh"
#include "TSPhysicsList.hh"
#include "TSActionInitialization.hh"

#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void message(RunManager* runmanager)
{
#ifdef G4MULTITHREADED
    runmanager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
    G4cout << "\n\n\t--> Running in multithreaded mode with "
           << runmanager->GetNumberOfThreads()
           << " threads\n\n" << G4endl;
#else
    // get rid of unused variable warning
    runmanager->SetVerboseLevel(runmanager->GetVerboseLevel());
    G4cout << "\n\n\t--> Running in serial mode\n\n" << G4endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{

    // Detect interactive mode (if no arguments) and define UI session
    //
    G4UIExecutive* ui = 0;
    if(argc == 1)
        ui = new G4UIExecutive(argc, argv);

    // Set the random seed
    CLHEP::HepRandom::setTheSeed(1245214UL);

    RunManager* runmanager = new RunManager();

    message(runmanager);

    runmanager->SetUserInitialization(new TSDetectorConstruction);

    runmanager->SetUserInitialization(new TSPhysicsList);

    runmanager->SetUserInitialization(new TSActionInitialization);

    runmanager->Initialize();


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
    if (!ui)
    {
        // batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[argc-1];
        UImanager->ApplyCommand(command+fileName);
    } else
    {
        // interactive mode
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
        delete ui;
    }

    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !

    delete visManager;
    delete runmanager;

    return 0;
}
