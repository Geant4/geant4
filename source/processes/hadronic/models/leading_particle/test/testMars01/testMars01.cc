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
// $Id: testMars01.cc,v 1.1 2001-12-13 14:58:40 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - testMars01 
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Mars01DetectorConstruction.hh"
#include "Mars01RunAction.hh"
#include "Mars01EventAction.hh"
#include "Mars01PrimaryGeneratorAction.hh"
#include "Mars01PhysicsList.hh"
#include "Mars01SteppingAction.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

int main(int argc,char** argv) {
  // Construct the default run manager

  G4RunManager* runManager = new G4RunManager;
  
  // set mandatory initialization classes
  runManager->SetUserInitialization(new Mars01DetectorConstruction);
  runManager->SetUserInitialization(new Mars01PhysicsList);

  // User Action classes
  runManager->SetUserAction(new Mars01RunAction);
  runManager->SetUserAction(new Mars01PrimaryGeneratorAction);
  runManager->SetUserAction(new Mars01EventAction);
  //runManager->SetUserAction(new Mars01SteppingAction);

  if(argc==1) {
    // G4UIterminal is a (dumb) terminal.
    G4UIsession* session = new G4UIterminal;
    session->SessionStart();
    delete session;
  }
  else {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  
  // job termination
  delete runManager;
  return 0;

}



