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
// $Id: test10.cc,v 1.4 2001-07-11 10:09:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT4 - test10 
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------
#include "Tst10DetectorConstruction.hh"
#include "Tst10RunAction.hh"
#include "Tst10EventAction.hh"
#include "Tst10PrimaryGeneratorAction.hh"
#include "Tst10PhysicsList.hh"
#include "Tst10VisManager.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

int main(int argc,char** argv) {

  // Set the default random engine to RanecuEngine
  RanecuEngine defaultEngine;
  HepRandom::setTheEngine(&defaultEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  runManager->SetUserInitialization(new Tst10DetectorConstruction);
  runManager->SetUserInitialization(new Tst10PhysicsList);

  // UserAction classes
  runManager->SetUserAction(new Tst10RunAction);
  runManager->SetUserAction(new Tst10EventAction);
  runManager->SetUserAction(new Tst10PrimaryGeneratorAction);

#ifdef G4VIS_USE
  G4cout << "Visualization init\n";
  G4VisManager* visManager = new Tst10VisManager();
  visManager -> Initialize ();
#endif

  // User interactions
  // Define (G)UI for interactive mode
  if(argc==1)
  { 
    // G4UIterminal is a (dumb) terminal.
    G4UIsession * session = new G4UIterminal;
    session->SessionStart();
    delete session;
  }
  else
  // Batch mode
  { 
    G4UImanager * UI = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

  delete runManager;
  return 0;
}

