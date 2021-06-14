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
// -------------------------------------------------------------
//      GEANT4 FF_Neutron_HP
//
//  Command line options:
//      -i ARG      : run in batch mode from script file ARG
//      -o ARG      : write output to file ARG
//                    (defaults to FF_Neutron_HP.out)
//      -n ARG      : multithreading with ARG number of threads
//                    (only works if Geant4 was compiled with
//                    multithreading enables)
//
//  =============== Begin Documentation Comments ===============
//!
//! \file       FissionFragment.cc
//! \author     B. Wendt (brycen.linn.wendt@cern.ch)
//! \date       June 06, 2014
//!
//! \brief      Main program of the FissionFragment example
//!
//! \details    Application demonstrating the Fission Fragment model as used
//!                 within the neutron_hp model. It demostrates the capability
//!                 for fission product containment by the cladding in a water
//!                 moderated sub-critical assembly.
//!             It could also be further extended to calculate the effective
//!                 multiplication factor of the subcritical assembly for
//!                 various loading schemes.
//!
//  ================ End Documentation Comments ================
//
//  Modified:
//
//  05-08-20                                              ARibon
//  Replaced deprecated HP environmental variables with UI commands
//  23-06-14                                              BWendt
//  Added check for NeutronHP fission generator environment variable
//
// -------------------------------------------------------------

#include "globals.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "QGSP_BIC_HP.hh"
#include "Randomize.hh"

#include "FFDetectorConstruction.hh"
#include "FFActionInitialization.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4ParticleHPManager.hh"

// Entry point
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int main(int argc, char* argv[])
{
    int result;
    unsigned int numberOfThreads = 1;

    G4String scriptFileName = "";
    G4String outputFileName = "FF_Neutron_HP.out";
    G4UImanager* UIManager = NULL;

    // Activate production of fission fragments in neutronHP
    G4ParticleHPManager::GetInstance()->SetProduceFissionFragments( true );

    char Force[] = "G4FORCENUMBEROFTHREADS";
    if(std::getenv(Force) != NULL) {
       char doNotForce[]="G4FORCENUMBEROFTHREADS=1";
       putenv(doNotForce);
    }

    // Indicate the example is starting
    G4cout << "####   Starting: " << argv[0] << "    ####" << G4endl;

    //  Parse the command line arguments, if any
    for(int i = 1;
        i < argc;
        i += 2)
    {
        // Ensure that this is actually a command
        if(argv[i][0] != '-')
        {
            G4cerr << G4endl << "!!!!" << G4endl;
            G4cerr << "!!!! Error in argument " << i + 1 << G4endl;
            G4cerr << "!!!! A command-line option was expected, but \""
                   << argv[i] << "\" was found" << G4endl;
            G4cerr << "!!!! " << argv[0] << " will now terminate" << G4endl;
            G4cerr << "!!!!" << G4endl << G4endl;

            return EXIT_FAILURE;
        }

        // Ensure that the command-line option has an associated argument
        if(!(i + 1 < argc))
        {
            G4cerr << G4endl << "!!!!" << G4endl;
            G4cerr << "!!!! Error in argument " << i + 2 << G4endl;
            G4cerr << "!!!! An argument was expected, but \"" << argv[i + 1]
                   << "\" was found" << G4endl;
            G4cerr << "!!!! Ensure that a space is used to separate the "
                      "option and argument" << G4endl;
            G4cerr << "!!!! " << argv[0] << " will now terminate" << G4endl;
            G4cerr << "!!!!" << G4endl << G4endl;

            return EXIT_FAILURE;
        }

        switch(argv[i][1])
        {
        case 'i':
            scriptFileName = "/control/execute ";
            scriptFileName.append(argv[i + 1]);
            break;

        case 'o':
            outputFileName = argv[i + 1];
            break;

        case 'n':
            result = sscanf(argv[i + 1],
                            "%u",
                            &numberOfThreads);
            if(result != 1)
            {
                G4cerr << G4endl << "!!!!" << G4endl;
                G4cerr << "!!!! Error in argument " << i + 2 << G4endl;
                G4cerr << "!!!! An positive number was expected, but \""
                       << argv[i + 1] << "\" was found" << G4endl;
                G4cerr << "!!!! " << argv[0] << " will now terminate"
                       << G4endl;
                G4cerr << "!!!!" << G4endl << G4endl;

                return EXIT_FAILURE;
            }
            break;

        default:
            G4cout << G4endl << "!!!!" << G4endl;
            G4cout << "!!!! Warning for command " << i + 1 << G4endl;
            G4cout << "!!!! \"" << argv[i] << "\" is not a valid command"
                   << G4endl;
            G4cout << "!!!! " << argv[0] << " will ignore \"" << argv[i]
                   << "\" and \"" << argv[i + 1] << "\"" << G4endl;
            G4cout << "!!!!" << G4endl << G4endl;
        }
    }

    // Instantiate G4UIExecutive if interactive mode
    G4UIExecutive* ui = nullptr;
    if (scriptFileName.length() == 0) {
      ui = new G4UIExecutive(argc, argv);
    }

    // Set the Random engine
    // A seed of 62737819 produced a maximum number of 67 events on the
    // author's system before timing out the nightly test
    const G4long seed = 62737819;
#ifndef NDEBUG
    G4cout << "MT RNG Seed: " << seed << G4endl;
#endif // NDEBUG
    G4Random::setTheEngine(new CLHEP::MTwistEngine(seed));

    // Initialize the multithreaded run manager
    auto* runManager = G4RunManagerFactory::CreateRunManager();
    runManager->SetNumberOfThreads(numberOfThreads);
    G4cout << "    Threads requested:    " << numberOfThreads << G4endl;
    G4cout << "    Threads started:      " << runManager->GetNumberOfThreads() << G4endl;

    // Set mandatory initialization classes
    runManager->SetUserInitialization(new FFDetectorConstruction());
    runManager->SetUserInitialization(new QGSP_BIC_HP());
    runManager->SetUserInitialization(new FFActionInitialization());

    // Initialize the Geant4 kernel
    runManager->Initialize();

    // Initialize visualization
    G4VisManager* visManager = new G4VisExecutive();
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    UIManager = G4UImanager::GetUIpointer();

    if(!ui)
    {
        // Batch mode
        UIManager->ApplyCommand(scriptFileName);
    } else
    {
        // Interactive mode
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


