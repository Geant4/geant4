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
// $Id: remsim.cc,v 1.12 2005/12/07 14:41:36 guatelli Exp $
// GEANT4 tag $Name: geant4-08-00 $

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
#include "G4VisExecutive.hh"
#endif
#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif 

int main(int argc,char** argv)
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  //  HepRandom :: setTheSeed(0);
  // Set mandatory initialization classes
  
  // Geometry
  RemSimDetectorConstruction* detector = new RemSimDetectorConstruction();
  runManager -> SetUserInitialization(detector);
  
  // Physics
  runManager->SetUserInitialization(new RemSimPhysicsList);

  // Set mandatory user action class

  // Primary particles
  RemSimPrimaryGeneratorAction* primary = new RemSimPrimaryGeneratorAction();
  runManager -> SetUserAction(primary);
 
  // Set optional user action class
  RemSimEventAction* event = new RemSimEventAction();
  runManager -> SetUserAction(event);

  RemSimRunAction* run = new RemSimRunAction();
  runManager -> SetUserAction(run);

  runManager -> SetUserAction(new RemSimSteppingAction(primary));

#ifdef G4VIS_USE
   // Visualisation
   G4VisManager* visManager = new G4VisExecutive;
   visManager -> Initialize();
#endif

#ifdef G4ANALYSIS_USE
 RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
 analysis -> SetFormat("hbook");
#endif
 
  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
   
  if(argc == 1)
    // Define (G)UI terminal for interactive mode  
    { 
      // G4UIterminal is a (dumb) terminal.
      G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif    

      UI -> ApplyCommand("/control/execute vis.mac");    
      session -> SessionStart();
      delete session;
    }
  else
    // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI -> ApplyCommand(command+fileName);
    }

#ifdef G4ANALYSIS_USE
 analysis -> finish();
#endif

#ifdef G4VIS_USE
 delete visManager;
#endif

  // job termination
  delete runManager;
 
  return 0;
}

