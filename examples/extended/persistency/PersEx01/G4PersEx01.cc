//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PersEx01.cc,v 1.1 2002-12-04 02:44:24 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN03 
//
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "Pers01VisManager.hh"
#endif

#include "Pers01DetectorConstruction.hh"
#include "Pers01PhysicsList.hh"
#include "Pers01PrimaryGeneratorAction.hh"
#include "Pers01RunAction.hh"
#include "Pers01EventAction.hh"
#include "Pers01SteppingAction.hh"
#include "Pers01SteppingVerbose.hh"

int main(int argc,char** argv) {

  // choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new Pers01SteppingVerbose);

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Pers01DetectorConstruction* detector = new Pers01DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Pers01PhysicsList);
  
 G4UIsession* session=0;
  
  if (argc==1)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_XM
      session = new G4UIXm(argc,argv);
#else           
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif
#endif
    }
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new Pers01VisManager;
  visManager->Initialize();
#endif
    
  // set user action classes
  runManager->SetUserAction(new Pers01PrimaryGeneratorAction(detector));
  runManager->SetUserAction(new Pers01RunAction);
  runManager->SetUserAction(new Pers01EventAction);
  runManager->SetUserAction(new Pers01SteppingAction);
  
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (session)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
      UI->ApplyCommand("/control/execute vis.mac");    
#ifdef G4UI_USE_XM
      // Customize the G4UIXm menubar with a macro file :
      UI->ApplyCommand("/control/execute visTutor/gui.mac");
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

