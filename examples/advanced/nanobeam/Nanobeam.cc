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
// -------------------------------------------------------------------
// $Id: Nanobeam.cc,v 1.8 2010-10-09 16:30:27 sincerti Exp $
// -------------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#ifdef G4UI_USE
  #include "G4UIExecutive.hh"
#endif

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"
#include "HistoManager.hh"

int main(int argc,char** argv) {

  // Choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
  
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysicsList);
  
  PrimaryGeneratorAction* primary = new PrimaryGeneratorAction(detector);
  runManager->SetUserAction(primary);
    
  HistoManager*  histo = new HistoManager();

  // Set user action classes
  RunAction* RunAct = new RunAction(detector,primary,histo);

  runManager->SetUserAction(RunAct);
  runManager->SetUserAction(new EventAction(RunAct));
  runManager->SetUserAction(new TrackingAction(RunAct)); 
  runManager->SetUserAction(new SteppingAction(RunAct,detector,primary,histo));
  
  // Initialize G4 kernel
  runManager->Initialize();
    
  // Get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  //
  remove ("nanobeam.root");
  
  if (argc==1)   // Define UI session for interactive mode.
  { 
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    UImanager->ApplyCommand("/control/execute default.mac");     
    ui->SessionStart();
    delete ui;
#endif
  }
     
  else           // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  delete runManager;

  return 0;
}
