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
/// \file exoticphysics/monopole/monopole.cc
/// \brief Main program of the exoticphysics/monopole example
//
// $Id: monopole.cc 66817 2013-01-12 16:16:08Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include "G4MonopolePhysics.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  //choose the Random engine
  //CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  //CLHEP::HepRandom::setTheEngine(new CLHEP::Ranlux64Engine);
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanluxEngine);

  //Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  //create physicsList
  // Physics List is defined via environment variable PHYSLIST
  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = factory.ReferencePhysList(); 
 
  // monopole physics is added
  G4MonopolePhysics * theMonopole = new G4MonopolePhysics();
  phys->RegisterPhysics(theMonopole);
    
  // visualization manager
#ifdef G4VIS_USE
  G4VisManager* visManager = 0;
#endif
        
  //get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Setup monopole
  G4String s = ""; 
  if(argc > 2) { s = argv[2]; }
  UImanager->ApplyCommand("/control/verbose 1");
  UImanager->ApplyCommand("/monopole/setup "+s);

  // mandator user classes
  DetectorConstruction* det = new DetectorConstruction();

  runManager->SetUserInitialization(det);  
  runManager->SetUserInitialization(phys);

  PrimaryGeneratorAction* kin = new PrimaryGeneratorAction(det);
  runManager->SetUserAction(kin);

  //user action classes
  RunAction* run;

  runManager->SetUserAction(run = new RunAction(det, kin));
  runManager->SetUserAction(new EventAction);
  runManager->SetUserAction(new TrackingAction(run));
  runManager->SetUserAction(new SteppingAction(run));

  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);    
    }
  else
    {  // interactive mode : define UI session
#ifdef G4VIS_USE
      //visualization manager
      visManager = new G4VisExecutive;
      visManager->Initialize();
#endif
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      ui->SessionStart();
      delete ui;
#endif
    }

  //job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
 
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
