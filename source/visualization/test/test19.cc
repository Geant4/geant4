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
// $Id: test19.cc,v 1.34 2010-10-26 16:17:44 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

// Usage: test19 [<session>] [<verbosity>]
// Without verbosity, verbosity=warnings.

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#include "globals.hh"
//#include "SharedSolidDetectorConstruction.hh"
#include "test19DetectorConstruction.hh"
#include "QGSP_BERT.hh"
#include "MyRunAction.hh"
#include "MyPrimaryGeneratorAction.hh"
#include "MyEventAction.hh"
#include "MySteppingAction.hh"

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4RunManager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
// G4VisExecutive is a G4VisManager that implements graphics system
// registration in the user domain.
#include "G4XXX.hh"
#include "G4XXXFile.hh"
#include "G4XXXStored.hh"
#ifdef G4VIS_USE_XXXSG
#include "G4XXXSG.hh"
#endif
#endif


#ifdef G4UI_USE_WIN32
#include <windows.h>
int WINAPI WinMain (
HINSTANCE hInstance,HINSTANCE hPrevInstance,LPSTR lpszCmdLine,int nCmdShow) {
#else
int main (int argc, char** argv) {
#endif

  G4String verbosityString("warnings");
#ifndef G4UI_USE_WIN32
  if ((argc >= 3)) verbosityString = argv[2];
#endif

#ifdef G4UI_USE
  // Choose (G)UI.
  G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#endif

  // Run manager
  G4cout << "RunManager is constructing...." << G4endl;
  G4RunManager * runManager = new G4RunManager;

  // User initialization classes
  runManager -> SetUserInitialization (new test19DetectorConstruction);
  //runManager -> SetUserInitialization (new SharedSolidDetectorConstruction);
  runManager -> SetUserInitialization (new QGSP_BERT);

  // UserAction classes.
  runManager -> SetUserAction (new MyRunAction);
  runManager -> SetUserAction (new MyPrimaryGeneratorAction);
  runManager -> SetUserAction (new MyEventAction);
  runManager -> SetUserAction (new MySteppingAction);

  //Initialize G4 kernel
  //runManager->Initialize();  // Do this with /run/initialize so that
			       // you can, optionally, choose detector
			       // (/test19det/detector N) first.

#ifdef G4VIS_USE
  // Instantiate and initialise Visualization Manager.
  G4VisManager* visManager = new G4VisExecutive;
  visManager -> SetVerboseLevel (verbosityString);
  visManager -> RegisterGraphicsSystem(new G4XXX);
  visManager -> RegisterGraphicsSystem(new G4XXXFile);
  visManager -> RegisterGraphicsSystem(new G4XXXStored);
#ifdef G4VIS_USE_XXXSG
  visManager -> RegisterGraphicsSystem(new G4XXXSG);
#endif
  visManager -> Initialize ();
#endif

  G4UImanager* UImanager = G4UImanager::GetUIpointer ();
#ifdef G4UI_USE_WIN32
  G4cout << "Reading win32.g4m file...." << G4endl;
  UImanager -> ApplyCommand ("/control/execute win32.g4m");
#else
  G4cout << "Reading test19.g4m file...." << G4endl;
  UImanager -> ApplyCommand ("/control/execute test19.g4m");
#endif

  G4cout << 
    "Choose a detector with /test19det/detector (or let default be"
    " constructed)."
       << G4endl;

#ifdef G4UI_USE
  // Start an interactive session.
  ui -> SessionStart();
#endif

#ifdef G4VIS_USE
  G4cout << "vis_test19: Deleting vis manager..." << G4endl;
  delete visManager;
  G4cout << "vis_test19: Vis manager deleted." << G4endl;
#endif
  G4cout << "vis_test19: Deleting run manager..." << G4endl;
  delete runManager;
  G4cout << "vis_test19: Run manager deleted." << G4endl;
  G4cout << "vis_test19: Deleting session..." << G4endl;
#ifdef G4UI_USE
  delete ui;
  G4cout << "vis_test19: Session deleted." << G4endl;
#endif

  return 0;
}
