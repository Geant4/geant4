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
// Code developed by:
//  S.Larsson
//
//    *********************
//    *                   *
//    *    PurgMag.cc     *
//    *                   *
//    *********************
//
// $Id: PurgMag.cc,v 1.2 2004/06/18 09:17:42 gunter Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// Comments: Main program for the Purgin Magnet example. 
//

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "PurgMagVisManager.hh"
#endif

#include "PurgMagDetectorConstruction.hh"
#include "PurgMagPhysicsList.hh"
#include "PurgMagPrimaryGeneratorAction.hh"
#include "PurgMagRunAction.hh"
#include "PurgMagEventAction.hh"
#include "PurgMagTrackingAction.hh"
#include "PurgMagSteppingAction.hh"
#include "PurgMagSteppingVerbose.hh"

int main(int argc,char** argv) {

  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  //PurgMag Verbose output class
  G4VSteppingVerbose::SetInstance(new PurgMagSteppingVerbose);
  
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  PurgMagDetectorConstruction* detector = new PurgMagDetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PurgMagPhysicsList);
  
  runManager->SetUserAction(new PurgMagPrimaryGeneratorAction(detector));
    
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new PurgMagVisManager;
  visManager->Initialize();
#endif

// output environment variables:
#ifdef G4ANALYSIS_USE
   G4cout << G4endl << G4endl << G4endl 
	  << " User Environment " << G4endl
	  << " Using AIDA 3.0 analysis " << G4endl;
# else
   G4cout << G4endl << G4endl << G4endl 
	  << " User Environment " << G4endl
	  << " G4ANALYSIS_USE environment variable not set, NO ANALYSIS " 
	  << G4endl;
#endif   

  // set user action classes
  PurgMagRunAction* RunAct = new PurgMagRunAction(detector);
  runManager->SetUserAction(RunAct);
  runManager->SetUserAction(new PurgMagEventAction(RunAct));
  runManager->SetUserAction(new PurgMagTrackingAction(RunAct)); 
  runManager->SetUserAction(new PurgMagSteppingAction(RunAct,detector));
  
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  


  if (argc==1)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
      G4UIsession * session = new G4UIterminal;
      UI->ApplyCommand("/control/execute vis.mac");    
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }


  // For bg jobs with 100000 events.
  /*
   int numberOfEvent = 100000;
   runManager->BeamOn(numberOfEvent);
  */


  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

