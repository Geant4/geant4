// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test201.cc,v 1.6 2000-06-08 12:49:20 barrand Exp $
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
G4std::ostream& g4cout = G4cout;
G4std::ostream& g4cerr = G4cerr;
#include "test201DetectorConstruction.hh"
#include "MyPhysicsList.hh"
#include "MyRunAction.hh"
#include "MyPrimaryGeneratorAction.hh"
#include "MyEventAction.hh"
#include "MySteppingAction.hh"

#include "G4StateManager.hh"

#ifdef G4UI_USE_TERMINAL
  #include "G4UIterminal.hh"
#endif
#ifdef G4UI_USE_GAG
  #include "G4UIGAG.hh"
#endif
#ifdef G4UI_USE_WO
  #include "G4UIWo.hh"
#endif
#ifdef G4UI_USE_XM
  #include "G4UIXm.hh"
#endif
#ifdef G4UI_USE_XAW
  #include "G4UIXaw.hh"
#endif
#ifdef G4UI_USE_WIN32
  #include "G4UIWin32.hh"
#endif

#include "G4RunManager.hh"

#ifdef G4VIS_USE
#include "MyVisManager.cc"
#endif

#ifdef G4UI_USE_WIN32
#include <windows.h>
int WINAPI WinMain (
HINSTANCE hInstance,HINSTANCE hPrevInstance,LPSTR lpszCmdLine,int nCmdShow) {
#else
int main (int argc, char** argv) {
#endif

  G4int Verbose = 0;
#ifndef G4UI_USE_WIN32
  if ((argc >= 3)) Verbose = atoi (argv[2]);
#endif

  // Choose (G)UI.
  G4UIsession* session;
#ifdef G4UI_USE_WIN32
  session = new G4UIWin32 (hInstance,hPrevInstance,lpszCmdLine,nCmdShow);
#else
  if (argc >= 2) {
    if (strcmp (argv[1], "dumb")==0)     session = new G4UIterminal;
#ifdef G4UI_USE_WO
    else if (strcmp (argv[1], "Wo")==0)  session = new G4UIWo (argc, argv);
#endif
#ifdef G4UI_USE_XM
    else if (strcmp (argv[1], "Xm")==0)  session = new G4UIXm (argc, argv);
#endif
#ifdef G4UI_USE_XAW
    else if (strcmp (argv[1], "Xaw")==0) session = new G4UIXaw (argc, argv);
#endif
#ifdef G4UI_USE_GAG
    else if (strcmp (argv[1], "gag")==0) session = new G4UIGAG ;
#endif
#ifdef G4UI_USE_GAG
    else                                 session = new G4UIGAG;
#else
    else                                 session = new G4UIterminal;
#endif
  } else {
#ifdef G4UI_USE_GAG
    session = new G4UIGAG;
#else
    session = new G4UIterminal;
#endif
  }
#endif
  G4UImanager::GetUIpointer()->SetSession(session);  //So that Pause works..

  // Run manager
  g4cout << "RunManager is constructing...." << G4endl;
  G4RunManager * runManager = new G4RunManager;

  // User initialization classes
  runManager -> SetUserInitialization (new test201DetectorConstruction);
  runManager -> SetUserInitialization (new MyPhysicsList);

  // UserAction classes.
  runManager -> SetUserAction (new MyRunAction);
  runManager -> SetUserAction (new MyPrimaryGeneratorAction);
  runManager -> SetUserAction (new MyEventAction);
  runManager -> SetUserAction (new MySteppingAction);

#ifdef G4VIS_USE
  G4VisManager* visManager = new MyVisManager;
  visManager -> SetVerboseLevel (Verbose);
  visManager -> Initialize ();
#endif

  G4UImanager* UI = G4UImanager::GetUIpointer ();

  g4cout << "Reading test201.g4m file...." << G4endl;
  UI -> ApplyCommand ("/control/execute test201.g4m");

  g4cout << 
    "Choose a detector with /test201/detector (or let default be"
    " constructed)."
       << G4endl;

  // Start an interactive session.
  session -> SessionStart();

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete session;
  delete runManager; // Should be last.

  return 0;
}

