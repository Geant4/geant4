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
// $Id: test59.cc,v 1.1 2009-11-27 16:06:28 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#define debug
#define visio

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "CHIPS.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"
#include "StackingAction.hh"
#include "HistoManager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

int main(int argc,char** argv)
{
  // Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  // Verbose output
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // Physics List
  G4VModularPhysicsList* physList = new CHIPS();
  runManager->SetUserInitialization( physList );

  // Initialization classes
  DetectorConstruction* detector;
  detector = new DetectorConstruction;      // ***** To be deleted
  runManager->SetUserInitialization(detector);

  HistoManager* histo = new HistoManager(); // ***** To be deleted

  // User actions ::::::::::::::::::::::::::::::::::::::::::::::::

  // Primary Generator
  PrimaryGeneratorAction* primary = new PrimaryGeneratorAction(detector);
  runManager->SetUserAction(primary);

  // Run Action
  RunAction* runact = new RunAction(detector,primary,histo);
  runManager->SetUserAction(runact);

  // Event Action
  EventAction* evact = new EventAction(runact, histo);
  runManager->SetUserAction(evact);

  // Track Action
  TrackingAction* trackingaction = new TrackingAction(detector, runact, evact, histo);
  runManager->SetUserAction(trackingaction);

  // Step Action
  SteppingAction* steppingaction = new SteppingAction(detector, runact, evact, histo);
  runManager->SetUserAction(steppingaction);
  
  // Stack Action
  StackingAction* stackingaction = new StackingAction(runact, evact, histo);
  runManager->SetUserAction(stackingaction);      
   
  // User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();  
 
  if (argc!=1)   // "batch" mode  
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }
  else           // define visualization and UI terminal for interactive mode
  { 
#ifdef visio
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif     
    G4UIsession * session = 0;
    session = new G4UIterminal(new G4UItcsh);      
    session = new G4UIterminal(); // for bush
    session->SessionStart();
    delete session;
#ifdef visio
    delete visManager;
#endif     
  }
    
  delete histo;
  delete runManager;

  return 0;
}
