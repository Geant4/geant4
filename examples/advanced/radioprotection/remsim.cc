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
// $Id: remsim.cc,v 1.3 2004-03-12 09:17:45 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN01 
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "RemSimDetectorConstruction.hh"
#include "RemSimPhysicsList.hh"
#include "RemSimPrimaryGeneratorAction.hh"
#include "RemSimEventAction.hh"
#include "RemSimRunAction.hh"
#include "RemSimSteppingAction.hh"
#ifdef G4VIS_USE
#include "RemSimVisManager.hh"
#endif
#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif 

#include "QGSP.hh"
#include "QGSC.hh" 
int main(int argc,char** argv)
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  RemSimDetectorConstruction* detector = new RemSimDetectorConstruction();
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new RemSimPhysicsList);

  // set mandatory user action class
  RemSimPrimaryGeneratorAction* primary = new RemSimPrimaryGeneratorAction();
  runManager->SetUserAction(primary);
 
  RemSimEventAction* event = new RemSimEventAction();
  runManager->SetUserAction(event);

  RemSimRunAction* run = new RemSimRunAction();
  runManager->SetUserAction(run);

  runManager->SetUserAction(new RemSimSteppingAction(primary,event,detector));

#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new RemSimVisManager;
  visManager->Initialize();
#endif
 
#ifdef G4ANALYSIS_USE
  RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
  analysis->book();
#endif   

  // Initialize G4 kernel
  // runManager->Initialize();

  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
   
if(argc==1)
  // Define (G)UI terminal for interactive mode  
  { 
    // G4UIterminal is a (dumb) terminal.
    G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif    

    UI->ApplyCommand("/control/execute vis1.mac");    
    session->SessionStart();
    delete session;
  }
  else
  // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }
 

#ifdef G4ANALYSIS_USE
//  RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
  analysis->finish();
#endif   

#ifdef G4VIS_USE
  delete visManager;
#endif

  // job termination
  delete runManager;
  return 0;
}


