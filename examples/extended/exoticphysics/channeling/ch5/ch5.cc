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
// gpaterno, October 2025
//
/// \file ch5.cc
/// \brief Main program of the ch5 example
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#include "G4ScoringManager.hh"
#include "G4AnalysisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "FTFP_BERT.hh"
#include "G4FastSimulationPhysics.hh"
#include "G4CoherentPairProductionPhysics.hh"

#include "Randomize.hh"
#include <ctime>
#include "G4Timer.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
    // Get current time
    G4Timer* theTimer = new G4Timer();
    theTimer->Start();
     
    //Set random number generator
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    G4String option_file = "random.in";
    std::ifstream fin(option_file);
    long random_seed = 0;
    if (fin.is_open()) {
        fin >> random_seed;
        fin.close();
    }
    random_seed += time(NULL);
    G4cout << "Random seed: " << random_seed << G4endl;
    CLHEP::HepRandom::setTheSeed(random_seed);
 
    //Use G4SteppingVerboseWithUnits
    G4int precision = 4;
    G4SteppingVerbose::UseBestUnit(precision);

    //Construct the run manager
    int vNumberOfThreads = 1;
    if (argc > 2) {
        vNumberOfThreads = atoi(argv[2]);
    }
    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    runManager->SetNumberOfThreads(vNumberOfThreads);
    G4cout << "### Using " << vNumberOfThreads << " threads ###" << G4endl;
                   
    //Activate UI-command based scorer
    G4ScoringManager* scManager = G4ScoringManager::GetScoringManager();
    scManager->SetVerboseLevel(0);

           
    //Set mandatory initialization classes    
    //Set the Geometry
    runManager->SetUserInitialization(new DetectorConstruction);

    //Physics list
    G4VModularPhysicsList* physicsList = new FTFP_BERT;
    // -- Create helper tool used to activate the fast simulation
    G4FastSimulationPhysics* fastSimulationPhysics = new G4FastSimulationPhysics();
    fastSimulationPhysics->BeVerbose();
    // -- activation of fast simulation for particles having fast simulation models
    // -- attached in the mass geometry
    fastSimulationPhysics->ActivateFastSimulation("e-");
    fastSimulationPhysics->ActivateFastSimulation("e+");
    fastSimulationPhysics->ActivateFastSimulation("pi-");
    fastSimulationPhysics->ActivateFastSimulation("pi+");
    fastSimulationPhysics->ActivateFastSimulation("mu-");
    fastSimulationPhysics->ActivateFastSimulation("mu+");
    fastSimulationPhysics->ActivateFastSimulation("proton");
    fastSimulationPhysics->ActivateFastSimulation("anti_proton");
    fastSimulationPhysics->ActivateFastSimulation("GenericIon");
    // -- Attach the fast simulation physics constructor to the physics list
    physicsList->RegisterPhysics(fastSimulationPhysics);
    /*
    //Coherent pair production model
    G4CoherentPairProductionPhysics* coherentPairProductionPhysics =
        new G4CoherentPairProductionPhysics();
    physicsList->RegisterPhysics(coherentPairProductionPhysics);
    */
    physicsList->SetVerboseLevel(1);
    runManager->SetUserInitialization(physicsList);

    //Set user action classes
    runManager->SetUserInitialization(new ActionInitialization());


    //Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
       
    if (argc != 1) {
        //Batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    } else {
        //Visualization manager
        G4VisManager* visManager = new G4VisExecutive;
        visManager->Initialize();

        //Define UI session for interactive mode        
        G4UIExecutive* ui = new G4UIExecutive(argc,argv);
        UImanager->ApplyCommand("/control/execute macros/init_vis.mac");
        if (ui->IsGUI()) UImanager->ApplyCommand("/control/execute macros/gui.mac");
        ui->SessionStart();
        delete ui;

        delete visManager;
    }
    
    //Job termination
    delete runManager;

    theTimer->Stop();
    G4cout << "Execution terminated" << G4endl;
    G4cout << (*theTimer) << G4endl;
    delete theTimer;
    
    return 0;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
