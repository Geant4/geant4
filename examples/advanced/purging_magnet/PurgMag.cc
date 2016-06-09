//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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
// $Id: PurgMag.cc,v 1.5 2006/06/29 16:05:47 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// Comments: Main program for the Purgin Magnet example. 
//

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
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
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
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
  G4VisManager* visManager = new G4VisExecutive;
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

