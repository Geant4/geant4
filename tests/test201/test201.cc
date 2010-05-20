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
// $Id: test201.cc,v 1.18 2010-05-20 18:10:48 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

// Usage: test201 dumb 5
//   Arguments have defaults: GAG 0
//   Second argument is a verbosity-controlling integer flag.

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#include "globals.hh"
std::ostream& g4cout = G4cout;
std::ostream& g4cerr = G4cerr;
#include "test201DetectorConstruction.hh"
#include "QGSP_BERT.hh"
#include "MyRunAction.hh"
#include "MyPrimaryGeneratorAction.hh"
#include "MyEventAction.hh"
#include "MySteppingAction.hh"
#include "G4UIsession.hh"
#include "G4UImanager.hh"
#include "G4StateManager.hh"
#include "G4RunManager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main (int argc, char** argv) {

  G4int Verbose = 0;
#ifndef G4UI_USE_WIN32
  if ((argc >= 3)) Verbose = atoi (argv[2]);
#endif

  // Run manager
  g4cout << "RunManager is constructing...." << G4endl;
  G4RunManager * runManager = new G4RunManager;

  // User initialization classes
  runManager -> SetUserInitialization (new test201DetectorConstruction);
  runManager -> SetUserInitialization (new QGSP_BERT);

  // UserAction classes.
  runManager -> SetUserAction (new MyRunAction);
  runManager -> SetUserAction (new MyPrimaryGeneratorAction);
  runManager -> SetUserAction (new MyEventAction);
  runManager -> SetUserAction (new MySteppingAction);

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager -> SetVerboseLevel (Verbose);
  visManager -> Initialize ();
#endif

  G4UImanager* UImanager = G4UImanager::GetUIpointer ();

  g4cout << "Reading test201.g4m file...." << G4endl;
  UImanager -> ApplyCommand ("/control/execute test201.g4m");

  g4cout << 
    "Choose a detector with /test201/detector (or let default be"
    " constructed)."
       << G4endl;

#ifdef G4UI_USE
  // Start an interactive session.
  G4UIExecutive* ui = new G4UIExecutive(argc, argv);
  ui -> SessionStart();
#endif

#ifdef G4UI_USE
  delete ui;
#endif
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager; // Should be last.

  return 0;
}
