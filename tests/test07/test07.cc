// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test07.cc,v 1.1 1999-01-08 16:35:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN02 
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
#include "G4UIterminal.hh"

#ifdef G4VIS_USE
#include "T07VisManager.hh"
#endif

#include "T07DetectorConstruction.hh"
#include "T07PhysicsList.hh"
#include "T07PrimaryGeneratorAction.hh"
#include "T07RunAction.hh"
#include "T07EventAction.hh"
#include "T07SteppingAction.hh"

int main(int argc,char** argv) {

#ifdef MAKEHBOOK

  theHbookManager.SetFilename("test07.hbook");

#endif

  // Set the default random engine to RanecuEngine
  RanecuEngine defaultEngine;
  HepRandom::setTheEngine(&defaultEngine);

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  T07DetectorConstruction* detector = new T07DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new T07PhysicsList);
  
#if defined(G4VIS_USE) && !defined(MAKEHBOOK)
  // visualization manager
  G4VisManager* visManager = new T07VisManager;
  visManager->Initialize();
#endif
    
  // set user action classes
  runManager->SetUserAction(new T07PrimaryGeneratorAction(detector));
  runManager->SetUserAction(new T07RunAction);
  runManager->SetUserAction(new T07EventAction);
  runManager->SetUserAction(new T07SteppingAction);
  
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc==1)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
      G4UIsession * session = new G4UIterminal;
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  // job termination
#if defined(G4VIS_USE) && !defined(MAKEHBOOK)
  delete visManager;
#endif
  delete runManager;

  return 0;
}




