// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: field01.cc,v 1.1 2001-03-27 16:20:06 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - TestF01 
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
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "F01VisManager.hh"
#endif

#include "F01DetectorConstruction.hh"
#include "F01ElectroMagneticField.hh"
#include "F01PhysicsList.hh"
#include "F01PrimaryGeneratorAction.hh"
#include "F01RunAction.hh"
#include "F01EventAction.hh"
#include "F01SteppingAction.hh"
#include "F01SteppingVerbose.hh"

int main(int argc,char** argv) 
{

  //choose the Random engine

  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class

  G4VSteppingVerbose::SetInstance(new F01SteppingVerbose);
  
  F01ElectroMagneticField* field = new F01ElectroMagneticField() ;
    
  // Construct the default run manager

  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes

  F01DetectorConstruction* detector;
  detector = new F01DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new F01PhysicsList(detector));
  
#ifdef G4VIS_USE

  // visualization manager

  G4VisManager* visManager = new F01VisManager;
  visManager->Initialize();

#endif 
 
  // set user action classes

  runManager->SetUserAction(new F01PrimaryGeneratorAction(detector));

  F01RunAction* runAction = new F01RunAction;

  runManager->SetUserAction(runAction);

  F01EventAction* eventAction = new F01EventAction(runAction);

  runManager->SetUserAction(eventAction);

  F01SteppingAction* steppingAction = new F01SteppingAction(detector,
                                                            eventAction, 
                                                            runAction);
  runManager->SetUserAction(steppingAction);
  
  //Initialize G4 kernel, physics tables ...

  runManager->Initialize();
    
  // get the pointer to the User Interface manager 

  G4UImanager* UI = G4UImanager::GetUIpointer();  
 
  if (argc==1)   // Define UI terminal for interactive mode  
  { 
     G4UIsession * session = new G4UIterminal;
     UI->ApplyCommand("/control/execute init.mac");    
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

