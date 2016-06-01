// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: myDetector.cc,v 1.1 1998/11/12 10:50:01 yhajime Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - GGEVis.cc
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIGAG.hh"
#include "MyDetectorConstruction.hh"
#include "MyPhysicsList.hh"
#include "MyPrimaryGeneratorAction.hh"

#ifdef G4VIS_USE
#include "MyVisManager.hh"
#endif


#include "g4templates.hh"

int main()
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  runManager->SetUserInitialization(new MyDetectorConstruction);
  runManager->SetUserInitialization(new MyPhysicsList);

#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new MyVisManager;
  visManager->Initialize();
#endif

  // set mandatory user action class
  runManager->SetUserAction(new MyPrimaryGeneratorAction);

  // Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

      G4UIsession * session = new G4UIGAG;
      session->SessionStart();
      delete session;
 

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;
  return 0;
}








