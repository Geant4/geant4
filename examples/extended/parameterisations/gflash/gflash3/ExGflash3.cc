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
/// \file ExGflash3.cc
/// \brief Main program of the parameterisations/gflash/gflash3 example
//
// Created by Joanna Weng 26.11.2004

// G4 includes 
#include "G4Types.hh"
#include "G4ios.hh"
#include "G4Timer.hh"
#include "G4UImanager.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

// my project 
#include "ExGflash3DetectorConstruction.hh"
#include "ExGflash3ParallelWorld.hh"
#include "ExGflashActionInitialization.hh"

#include "FTFP_BERT.hh"
#include "G4FastSimulationPhysics.hh"
#include "G4ParallelWorldPhysics.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{   
  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // timer to see GFlash performance
  G4Timer timer;
  timer.Start();
  
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
  G4cout<<"|                                                       |"<<G4endl;
  G4cout<<"|          This is an example of Shower                 |"<<G4endl;
  G4cout<<"|          Parameterization with GFLASH                 |"<<G4endl;
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
  


#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(1);
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
  G4cout<<"|              Constructing MT run manager              |"<<G4endl;
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
#else
  G4RunManager * runManager = new G4RunManager;
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
  G4cout<<"|        Constructing sequential run manager            |"<<G4endl;
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
#endif
  
  // UserInitialization classes (mandatory)
  G4cout<<"# GFlash Example: Detector Construction"<<G4endl;    
  auto detector = new ExGflash3DetectorConstruction();
  detector->RegisterParallelWorld(new ExGflash3ParallelWorld("parallelWorld"));
  runManager->SetUserInitialization(detector);

  // G4cout<<"# GFlash Example: Physics "<<G4endl;
  // -- Select a physics list:
  G4VModularPhysicsList* physicsList = new FTFP_BERT();
  // -- Create a fast simulation physics constructor, used to augment
  // -- the above physics list to allow for fast simulation:
  G4FastSimulationPhysics* fastSimulationPhysics = new G4FastSimulationPhysics();
  // -- We now configure the fastSimulationPhysics object.
  // -- The gflash model (GFlashShowerModel, see ExGflashDetectorConstruction.cc)
  // -- is applicable to e+ and e- : we augment the physics list for these
  // -- particles (by adding a G4FastSimulationManagerProcess with below's
  // -- calls), this will make the fast simulation to be activated:
  fastSimulationPhysics->ActivateFastSimulation("e-");
  fastSimulationPhysics->ActivateFastSimulation("e+");
  // -- Register this fastSimulationPhysics to the physicsList,
  // -- when the physics list will be called by the run manager
  // -- (will happen at initialization of the run manager)
  // -- for physics process construction, the fast simulation
  // -- configuration will be applied as well.
  physicsList->RegisterPhysics( new G4ParallelWorldPhysics("parallelWorld") );
  physicsList->RegisterPhysics( fastSimulationPhysics );
  runManager->SetUserInitialization(physicsList);

  // Action initialization:
  runManager->SetUserInitialization(new ExGflashActionInitialization);

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/run/verbose 0");
  runManager->Initialize();
  UImanager->ApplyCommand("/Step/Verbose 0");
  
  if (ui)   // Define UI terminal for interactive mode
  { 
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  }
  else           // Batch mode
  { 
    G4String s=*(argv+1);
    UImanager->ApplyCommand("/control/execute "+s);
  }
  
  delete visManager;
  delete runManager;
  
  timer.Stop();
  G4cout << G4endl;
  G4cout << "******************************************";
  G4cout << G4endl;
  G4cout << "Total Real Elapsed Time is: "<< timer.GetRealElapsed();
  G4cout << G4endl;
  G4cout << "Total System Elapsed Time: " << timer.GetSystemElapsed();
  G4cout << G4endl;
  G4cout << "Total GetUserElapsed Time: " << timer.GetUserElapsed();
  G4cout << G4endl;
  G4cout << "******************************************";
  G4cout << G4endl;
  
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
