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
// $Id: test17.cc,v 1.8 2001-10-30 22:12:22 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - test17 
//
// --------------------------------------------------------------
// Comments
//     
// 18.08.2000 V.Ivanchenko clean up visualisation and dummy output
//   
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"

#include "Test17DetectorConstruction.hh"
#include "Test17PhysicsList.hh"
#include "Test17PrimaryGeneratorAction.hh"
#include "Test17RunAction.hh"
#include "Test17EventAction.hh"
#include "Test17SteppingAction.hh"

int main(int argc,char** argv) {

  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Test17DetectorConstruction* detector;
  detector = new Test17DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Test17PhysicsList(detector));
   
  // set user action classes
  runManager->SetUserAction(new Test17PrimaryGeneratorAction(detector));
  Test17RunAction* runaction = new Test17RunAction;
  runManager->SetUserAction(runaction);

  Test17EventAction* eventaction = new Test17EventAction(runaction);
  runManager->SetUserAction(eventaction);

  Test17SteppingAction* steppingaction = new Test17SteppingAction(detector,
                                               eventaction, runaction);
  runManager->SetUserAction(steppingaction);
  
  //Initialize G4 kernel
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
  delete runManager;

  return 0;
}

