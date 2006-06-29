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
// $Id: test10.cc,v 1.7 2006-06-29 21:38:10 gunter Exp $
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
  CLHEP::RanecuEngine defaultEngine;
  CLHEP::HepRandom::setTheEngine(&defaultEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  Tst10DetectorConstruction* detector =  new Tst10DetectorConstruction();
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Tst10PhysicsList);

  // UserAction classes
  runManager->SetUserAction(new Tst10RunAction);
  runManager->SetUserAction(new Tst10EventAction);
  runManager->SetUserAction(new Tst10PrimaryGeneratorAction(detector));

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

