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
//
// $Id: test32.cc,v 1.2 2006-06-29 21:58:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - testTst32 
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Tst32DetectorConstruction.hh"
#include "Tst32RunAction.hh"
#include "Tst32EventAction.hh"
#include "Tst32PrimaryGeneratorAction.hh"
#include "Tst32PhysicsList.hh"
#include "Tst32SteppingAction.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

int main(int argc,char** argv) {
  // Construct the default run manager

  G4RunManager* runManager = new G4RunManager;
  
  // set mandatory initialization classes
  runManager->SetUserInitialization(new Tst32DetectorConstruction);
  runManager->SetUserInitialization(new Tst32PhysicsList);

  // User Action classes
  runManager->SetUserAction(new Tst32RunAction);
  runManager->SetUserAction(new Tst32PrimaryGeneratorAction);
  runManager->SetUserAction(new Tst32EventAction);
  //runManager->SetUserAction(new Tst32SteppingAction);

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



