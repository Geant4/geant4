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
<<<<<<< HEAD:examples/extended/parameterisations/gflash/ExGflash.cc
// $Id: ExGflash.cc 94396 2015-11-13 13:37:16Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/extended/parameterisations/gflash/gflash1/ExGflash1.cc
//
/// \file ExGflash1.cc
/// \brief Main program of the parameterisations/gflash/gflash1 example
//
// Created by Joanna Weng 26.11.2004

// G4 includes
#include "G4Types.hh"
#include "G4ios.hh"
#include "G4Timer.hh"
#include "G4UImanager.hh"

#include "G4RunManagerFactory.hh"

// my project
#include "ExGflash1DetectorConstruction.hh"
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/extended/parameterisations/gflash/gflash1/ExGflash1.cc
#include "ExGflashActionInitialization.hh"

#include "FTFP_BERT.hh"
#include "G4VModularPhysicsList.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

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

  auto* runManager = G4RunManagerFactory::CreateRunManager();
  runManager->SetNumberOfThreads(1);

  // UserInitialization classes (mandatory)
  G4cout<<"# GFlash Example: Detector Construction"<<G4endl;
  runManager->SetUserInitialization(new ExGflash1DetectorConstruction);

  // G4cout<<"# GFlash Example: Physics "<<G4endl;
  G4VModularPhysicsList* physicsList = new FTFP_BERT();
  physicsList->RegisterPhysics(new ExGflashPhysics());
  runManager->SetUserInitialization(physicsList);

  // Action initialization:
  runManager->SetUserInitialization(new ExGflashActionInitialization);

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/run/verbose 0");
  runManager->Initialize();
  UImanager->ApplyCommand("/Step/Verbose 0");

  if (ui)   // Define UI terminal for interactive mode
  {
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
#endif
  }
  else           // Batch mode
  {
    G4String s=*(argv+1);
    UImanager->ApplyCommand("/control/execute "+s);
  }

  delete visManager;
  #endif  
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
