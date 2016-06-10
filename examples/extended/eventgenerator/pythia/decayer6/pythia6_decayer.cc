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
// $Id: pythia6_decayer.cc 72244 2013-07-12 08:49:56Z gcosmo $
//
/// \file eventgenerator/pythia/decayer6/pythia6_decayer.cc
/// \brief Main program of the pythia6Decayer example

#include "P6DExtDecayerPhysics.hh"
//#include "P6DPhysicsList.hh"

#include "ExG4DetectorConstruction01.hh"
#include "ExG4PrimaryGeneratorAction01.hh"
#include "ExG4RunAction01.hh"
#include "ExG4EventAction01.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ThreeVector.hh"
#include "QGSP_BERT.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv)
{
  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(new ExG4DetectorConstruction01);
  
  //
/*  
  G4VUserPhysicsList* physics = new P6DPhysicsList;
  runManager->SetUserInitialization(physics);
*/
  G4VModularPhysicsList* physicsList = new QGSP_BERT;
  physicsList->RegisterPhysics(new P6DExtDecayerPhysics());
  runManager->SetUserInitialization(physicsList);

    
  // Set user action classes
  //
  runManager->SetUserAction(
    new ExG4PrimaryGeneratorAction01("B-", 50.*MeV));
    // B- meson has not defined decay in Geant4
  runManager->SetUserAction(new ExG4RunAction01);
  runManager->SetUserAction(new ExG4EventAction01);
  //
  //G4UserSteppingAction* stepping_action =
  //                  new SteppingAction(detector, event_action);
  //runManager->SetUserAction(stepping_action);
  
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();      
  
  if (argc!=1) {
    // batch mode{
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);    
  }
  else {  // interactive mode : define UI session
#ifdef G4UI_USE
#ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac"); 
#else
    UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
    ui->SessionStart();
    delete ui;
#ifdef G4VIS_USE
    delete visManager;
#endif                
#endif
    }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  delete runManager;

  return 0;
}
