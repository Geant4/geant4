// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: lArCal.cc,v 1.2 2002-10-02 20:54:53 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN03 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "FCALVisManager.hh"
#endif

#include "FCALTestbeamSetup.hh"
#include "FCALPhysicsList.hh"
#include "FCALPrimaryGeneratorAction.hh"
#include "FCALRunAction.hh"
#include "FCALTBEventAction.hh"
#include "FCALSteppingAction.hh"
#include "FCALSteppingVerbose.hh"

int main(int argc,char** argv) {

  // choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new FCALSteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  FCALTestbeamSetup* detector = new FCALTestbeamSetup;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new FCALPhysicsList);
  
 G4UIsession* session=0;
  
  if (argc==1)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_XM
      session = new G4UIXm(argc,argv);
#else
      session = new G4UIterminal;
#endif
    }
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new FCALVisManager;
  visManager->Initialize();
#endif
    
  // set user action classes
  runManager->SetUserAction(new FCALPrimaryGeneratorAction());
  runManager->SetUserAction(new FCALRunAction);
  FCALSteppingAction* StepAction = new FCALSteppingAction;
  runManager->SetUserAction(new FCALTBEventAction(StepAction));
  runManager->SetUserAction(StepAction);
  
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (session)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
      UI->ApplyCommand("/control/execute prerunN03.mac");    
#ifdef G4UI_USE_XM
      // Customize the G4UIXm menubar with a macro file :
      UI->ApplyCommand("/control/execute gui.mac");
#endif
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
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

