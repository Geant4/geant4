// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// --------------------------------------------------------------
//      GEANT 4 - test03
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
//
// File:        test03.cc
// Description: Test of Continuous Process G4Cerenkov
//              -- Generation Cerenkov Photons --
//              -- Transport of optical Photons --
// Version:     5.0
// Created:     1996-04-30
// Author:      Juliet Armstrong
// Updated:     1998-06-03 by Peter Gumplinger for example Tst03
// mail:        gum@triumf.ca
// Updated:     1998-08-07 by John Allison - copied from example Tst03 to test03.
//
// --------------------------------------------------------------

#include "Tst03RunAction.hh"
#include "Tst03DetectorConstruction.hh"
#include "Tst03PrimaryGeneratorAction.hh"
#include "Tst03PhysicsList.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4ios.hh"
#include <stdlib.h>

int main(int argc,char** argv) {

  // Set the default random engine to RanecuEngine
  RanecuEngine defaultEngine;
  HepRandom::setTheEngine(&defaultEngine);

  // Seed the random number generator manually
  // -----------------------------------------

  G4long myseed = 345354;

  HepRandom::setTheSeed(myseed);

  // Run manager

  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes - mandatory

  runManager-> SetUserInitialization(new Tst03DetectorConstruction);
  runManager-> SetUserInitialization(new Tst03PhysicsList);

  // UserAction classes - optional

  runManager->SetUserAction(new Tst03RunAction);
  runManager->SetUserAction(new Tst03PrimaryGeneratorAction);

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

  return EXIT_SUCCESS;
}

void LoopUntilPressEnter()
{
        char ch;

        G4cout << "Press <Enter> to continue ... " << endl;
        while ( cin.get(ch) )
        {
                if (ch == '\n') break;
        }
        G4cout << endl;
}

