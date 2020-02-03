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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Ramos-Mendez, et. Al. Flagged uniform particle splitting for track  
//                  structure of proton and light ions 
//                  (accepted in Phys. Med. Biol) 
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file splitting.cc
/// \brief Main program of the splitting example

#include "G4Types.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4UImanager.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
    // Instantiate G4UIExecutive if interactive mode
    G4UIExecutive* ui = nullptr;
    if ( argc == 1 ) {
      ui = new G4UIExecutive(argc, argv);
    }

#ifdef G4MULTITHREADED
    G4MTRunManager* runManager= new G4MTRunManager;
    G4int nThreads = 2;
    runManager->SetNumberOfThreads(nThreads);
#else
    G4RunManager* runManager = new G4RunManager();
#endif
  
    // Set user classes
    //
    runManager->SetUserInitialization(new DetectorConstruction());
    PhysicsList* physList = new PhysicsList();
    runManager->SetUserInitialization(physList);
    runManager->SetUserInitialization(new ActionInitialization(physList));

    // Initialize G4 kernel
    G4UImanager* UI = G4UImanager::GetUIpointer();

    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    if (!ui)   // batch mode
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UI->ApplyCommand(command+fileName);
    }

    else           //define visualization and UI terminal for interactive mode
    {
        ui->SessionStart();
        delete ui;
    }

    // Job termination
    delete visManager;
    delete runManager;

    return 0;
}

