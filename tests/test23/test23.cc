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
// $Id: test23.cc,v 1.4 2004-03-18 11:02:24 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  Test of the G4QCaptureAtRest CHIPS process in GEANT4

#include "G4ios.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4RunManager.hh"

#include "Tst23DetectorConstruction.hh"
#include "Tst23RunAction.hh"
#include "Tst23PrimaryGeneratorAction.hh"
#include "Tst23PhysicsList.hh"
#include "Tst23SteppingAction.hh"
#include "Test23EventAction.hh"
#include "Test23TrackingAction.hh"


int main(int argc,char** argv) {

  // Set the default random engine to RanecuEngine
  RanecuEngine defaultEngine;
  HepRandom::setTheEngine(&defaultEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  runManager->SetUserInitialization(new Tst23DetectorConstruction);
  runManager->SetUserInitialization(new Tst23PhysicsList);

  //set additional user action classes
  Test23EventAction* eventAction       = new Test23EventAction();
  Test23TrackingAction* trackingAction = new Test23TrackingAction();

  // UserAction classes
  runManager->SetUserAction(new Tst23RunAction); // (delete equivalent)
  runManager->SetUserAction(eventAction);        // (delete equivalent)
  runManager->SetUserAction(trackingAction);     // (delete equivalent)
  runManager->SetUserAction(new Tst23PrimaryGeneratorAction); // (delete equivalent)
  //runManager->SetUserAction(new Tst23SteppingAction); // (delete equivalent)

  if(argc==1)
  {
    // G4UIterminal is a (dumb) terminal.
    G4UIsession* session = new G4UIterminal;
    session->SessionStart();
    delete session;
  }
  else
  {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  trackingAction->PrintResult(); // Print the result
  G4cerr<<"Test23 is done"<<G4endl;

  delete runManager;
  //G4cerr<<"Test23:RunManager is deleted: if not reached, don't delete UserAct's"<<G4endl;
  return 0;
}

