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
/// \file channeling/channeling.cc
/// \brief Main program of the channeling example
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4ScoringManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "ExExChDetectorConstruction.hh"

#ifdef G4MULTITHREADED
#include "ExExChUserActionInitialization.hh"
#else
#include "ExExChPrimaryGeneratorAction.hh"
#include "ExExChTrackingAction.hh"
#include "ExExChStackingAction.hh"
#include "ExExChEventAction.hh"
#include "ExExChRunAction.hh"
#endif

#include "ExExChPhysicsList.hh"
#include "QGSP_BERT.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
    // Construct the default run manager
#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    if(argv[2]){
        if(atoi(argv[2])>0){
            runManager->SetNumberOfThreads(atoi(argv[2]));
        }
    }
    G4cout << "MT MODE ON " << runManager->GetNumberOfThreads() << G4endl;
#else
    G4RunManager* runManager = new G4RunManager;
    G4cout << "MT MODE OFF" << G4endl;
#endif
    
    // Activate UI-command base scorer
    G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
    scManager->SetVerboseLevel(0);
    
    // Choose the Random engine
#ifndef G4MULTITHREADED
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
#endif
    
    // Set mandatory initialization classes
    G4VUserDetectorConstruction* detector =
        new ExExChDetectorConstruction;
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new ExExChPhysicsList());

#ifndef G4MULTITHREADED
    // Set user action classes
    runManager->SetUserAction(new ExExChPrimaryGeneratorAction());
    runManager->SetUserAction(new ExExChEventAction());
    runManager->SetUserAction(new ExExChStackingAction());
    runManager->SetUserAction(new ExExChTrackingAction());
    runManager->SetUserAction(new ExExChRunAction());
#else
    runManager->SetUserInitialization(
                new ExExChUserActionInitialization());
#endif
    
    // Get the pointer to the User Interface manager
    G4UImanager* UI = G4UImanager::GetUIpointer();  
    
    if(argc!=1) {
        // Batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UI->ApplyCommand(command+fileName);
    }
    
    else {
        // Define visualization and UI terminal for interactive mode
#ifdef G4VIS_USE
        G4VisManager* visManager = new G4VisExecutive;
        visManager->Initialize();
#endif
        
#ifdef G4UI_USE
        G4UIExecutive * ui = new G4UIExecutive(argc,argv);
        ui->SessionStart();
        delete ui;
#endif
        
#ifdef G4VIS_USE
        delete visManager;
#endif
    }
    
    // Job termination
    delete runManager;
    
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
