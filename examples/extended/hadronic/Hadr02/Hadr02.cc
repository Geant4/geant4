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
/// \file hadronic/Hadr02/Hadr02.cc
/// \brief Main program of the hadronic/Hadr02 example
//
//
// -------------------------------------------------------------
//      GEANT4 Hadr02
//
//  Application demonstrating Geant4 hadronic cross sections
//
//  Author: V.Ivanchenko 20 June 2008
//
//  Modified: 
//
// -------------------------------------------------------------
//
//
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "StackingAction.hh"
#include "HistoManager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "UrQMD.hh"
#include "CRMC_FTFP_BERT.hh"


int main(int argc,char** argv) {

  //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc,argv);

  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());

  //Construct the default run manager
  G4RunManager * runManager = new G4RunManager();

  //set mandatory initialization classes
  runManager->SetUserInitialization(new DetectorConstruction());

  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = 0;

  // default Physics List for this example
  G4String physName = "QBBC";

  // Physics List name defined via 2nd argument
  if (argc==3) { physName = argv[2]; }
  else {
    char* path = std::getenv("PHYSLIST");
    if (path) { physName = G4String(path); }
  }
  if ( physName == "UrQMD" ) { 
    phys = new UrQMD; 
  } else if ( physName == "CRMC_FTFP_BERT" ) {
    phys = new CRMC_FTFP_BERT;
  } else { 
    phys = factory.GetReferencePhysList( physName ); 
  }

  // Physics List is defined via environment variable PHYSLIST
  if(!phys) {
    G4cout << "Hadr02 FATAL ERROR: Physics List is not defined"
           << G4endl;
    return 1;
  }
  runManager->SetUserInitialization(phys);
  HistoManager::GetPointer()->SetPhysicsList(phys);

  runManager->SetUserAction(new PrimaryGeneratorAction());

  //set user action classes
  runManager->SetUserAction(new RunAction());
  runManager->SetUserAction(new EventAction());
  runManager->SetUserAction(new StackingAction());

  //initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  //get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if (ui)  {
   //interactive mode
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

