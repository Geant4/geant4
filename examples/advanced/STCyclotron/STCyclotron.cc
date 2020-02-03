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
// --------------------------------------------------------------
//                 GEANT 4 - Solid Target Cyclotron example
// --------------------------------------------------------------
//
// Code developed currently by:
//  F. Poignant & S.Guatelli 
// Code also developed in the past by:
//  S. Penfold
//
// file STCyclotron.cc
//
#include "G4MTRunManager.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh" 
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "STCyclotronAnalysis.hh"

#include "STCyclotronActionInitialization.hh"
#include "STCyclotronDetectorConstruction.hh"
#include "STCyclotronPhysicsList.hh"
#include "STCyclotronPrimaryGeneratorAction.hh"
#include "STCyclotronRunAction.hh"

int main(int argc, char** argv)
{
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(4); // Is equal to 2 by default
#else
 G4RunManager* runManager = new G4RunManager;
#endif

  //Set mandatory initialization classes
  STCyclotronDetectorConstruction* det = new STCyclotronDetectorConstruction();
  runManager->SetUserInitialization(det);  
 
  STCyclotronPhysicsList* physList = new STCyclotronPhysicsList(det);
  runManager->SetUserInitialization(physList);

  // User action initialization
  STCyclotronActionInitialization* actions = new STCyclotronActionInitialization(det);
  runManager->SetUserInitialization(actions);

  //Set the random number seed
/* 
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  G4long seed=time(0);
  CLHEP::HepRandom::setTheSeed(seed);
  G4cout << "***********************" << G4endl;
  G4cout << "*** Seed: " << CLHEP::HepRandom::getTheSeed() << " ***" << G4endl;
  G4cout << "***********************" << G4endl;
  //getchar();
 */
 
  G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  //initialize : beam parameters + simulation parameters (geometry)
  UImanager->ApplyCommand("/control/execute Macro/init_parameters.mac");

  if (argc!=1) {
    // batch mode
    //to apply the command : sahmri_simulation run.mac
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  
  else {
    // interactive mode : define UI session

    G4UIExecutive* ui = new G4UIExecutive(argc, argv);

    UImanager->ApplyCommand("/control/execute Macro/Vis/init_vis.mac"); 
    //UImanager->ApplyCommand("/control/execute Macro/init.mac"); 

    if (ui->IsGUI())
      UImanager->ApplyCommand("/control/execute Macro/GUI/gui.mac");
     ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;

  return 0;

}
