// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test19.cc,v 1.8 2000-07-03 10:32:24 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

// Usage: test19 dumb 5
//   Arguments have defaults: GAG 0
//   Second argument is a verbosity-controlling integer flag.

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#include "globals.hh"
#include "test19DetectorConstruction.hh"
#include "MyPhysicsList.hh"
#include "MyRunAction.hh"
#include "MyPrimaryGeneratorAction.hh"
#include "MyEventAction.hh"
#include "MySteppingAction.hh"

#include "G4VTreeGraphicsScene.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4ModelingParameters.hh"

#include "G4StateManager.hh"

#ifdef G4UI_USE_TERMINAL
  #include "G4UIterminal.hh"
  #include "G4UItcsh.hh"
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

#include "MyVisManager.cc"

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
    if (strcmp (argv[1], "dumb")==0)     session =
					   new G4UIterminal(new G4UItcsh);
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
    else                                 session =
					   new G4UIterminal(new G4UItcsh);
  }
  else                                   session =
					   new G4UIterminal(new G4UItcsh);
#endif
  G4UImanager::GetUIpointer()->SetSession(session);  //So that Pause works..

  // Run manager
  G4cout << "RunManager is constructing...." << G4endl;
  G4RunManager * runManager = new G4RunManager;

  // User initialization classes
  runManager -> SetUserInitialization (new test19DetectorConstruction);
  runManager -> SetUserInitialization (new MyPhysicsList);

  // UserAction classes.
  runManager -> SetUserAction (new MyRunAction);
  runManager -> SetUserAction (new MyPrimaryGeneratorAction);
  runManager -> SetUserAction (new MyEventAction);
  runManager -> SetUserAction (new MySteppingAction);

  //Initialize G4 kernel
  runManager->Initialize();

  // Demonstration of G4VTreeGraphicsScene (DTREE).
  G4VTreeGraphicsScene treeDump;
  G4VPhysicalVolume* world =
    G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking()->GetWorldVolume();
  G4PhysicalVolumeModel worldModel(world);
  const G4ModelingParameters modelingParams;  // Default modeling parameters.
  worldModel.SetModelingParameters(&modelingParams);
  worldModel.DescribeYourselfTo(treeDump);

  // Instantiate and initialise Visualization Manager.
  G4VisManager* visManager = new MyVisManager;
  visManager -> SetVerboseLevel (Verbose);
  visManager -> Initialize ();

  G4UImanager* UI = G4UImanager::GetUIpointer ();

#ifdef G4UI_USE_WIN32
  G4cout << "Reading win32.g4m file...." << G4endl;
  UI -> ApplyCommand ("/control/execute win32.g4m");
#else
  G4cout << "Reading test19.g4m file...." << G4endl;
  UI -> ApplyCommand ("/control/execute test19.g4m");
#endif

  G4cout << 
    "Choose a detector with /test19det/detector (or let default be"
    " constructed)."
       << G4endl;

  // Start an interactive session.
  session -> SessionStart();

  G4cout << "vis_test19: Deleting vis manager..." << G4endl;
  delete visManager;
  G4cout << "vis_test19: Vis manager deleted." << G4endl;
  G4cout << "vis_test19: Deleting run manager..." << G4endl;
  delete runManager;
  G4cout << "vis_test19: Run manager deleted." << G4endl;
  G4cout << "vis_test19: Deleting session..." << G4endl;
  delete session;
  G4cout << "vis_test19: Session deleted." << G4endl;

  return 0;
}
