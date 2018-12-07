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
/// \file medical/electronScattering2/electronScattering2.cc
/// \brief Main program of the medical/electronScattering2 example
#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

#include "ElectronBenchmarkDetector.hh"
#include "PhysicsList.hh"
#include "ElectronActionInitialization.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

    //detect interactive mode (if no arguments) and define UI session
    G4UIExecutive* ui = nullptr;
    if (argc == 1) ui = new G4UIExecutive(argc,argv);

    // Parse the arguments
    G4String outputFile = "output.csv";
    G4String startingSeed = "1";
    G4String macroFile = "None";
    if (argc > 1) macroFile = argv[1];
    if (argc > 2) startingSeed = argv[2];
    if (argc > 3) outputFile = argv[3];
    G4cout << "Starting run with" << G4endl;
    G4cout << "Macro File    : " << macroFile << G4endl;
    G4cout << "Starting Seed : " << startingSeed << G4endl;
    G4cout << "Output File   : " << outputFile << G4endl;

    // Instantiate the run manager
#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
#else
    G4RunManager* runManager = new G4RunManager;
#endif

    // Instantiate the random engine
    G4Random::setTheEngine(new CLHEP::MTwistEngine);

    // Convert the starting seed to integer and feed it to the random engine
    unsigned startingSeedInt;
    std::istringstream is(startingSeed);
    is >> startingSeedInt;
    G4Random::setTheSeed(startingSeedInt);

    // Instantiate the geometry
    runManager->SetUserInitialization(new ElectronBenchmarkDetector);

    // Instantiate the physics list (in turn calls one of the choices of
    // physics list)
    runManager->SetUserInitialization(new PhysicsList);

    // set user action classes
    runManager->SetUserInitialization(
        new ElectronActionInitialization(outputFile));

    //initialize visualization
    G4VisManager* visManager = nullptr;

    //get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (ui)  {
     //interactive mode
     visManager = new G4VisExecutive;
     visManager->Initialize();
     ui->SessionStart();
     delete ui;
    }
    else  {
    //batch mode
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UImanager->ApplyCommand(command+fileName);
    }

    //job termination
    delete visManager;
    delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
