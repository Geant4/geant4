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
// $Id: testAssemblyVolume.cc,v 1.3 2001-07-11 09:59:22 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT4 - testAssemblyVolume 
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "TstVADetectorConstruction.hh"
#include "TstVARunAction.hh"
#include "TstVAEventAction.hh"
#include "TstVAPrimaryGeneratorAction.hh"
#include "TstVASteppingAction.hh"
#include "TstVAPhysicsList.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4RunManager.hh"
#ifdef G4VIS_USE
#include "TstVAVisManager.hh"
#endif

#include "G4ios.hh"

int main(int argc,char** argv) {

  // Set the default random engine to RanecuEngine
  RanecuEngine defaultEngine;
  HepRandom::setTheEngine(&defaultEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  runManager->SetUserInitialization(new TstVADetectorConstruction);
  runManager->SetUserInitialization(new TstVAPhysicsList);

  #ifdef G4VIS_USE
    // Visualization, if you choose to have it!
    G4VisManager* visManager = new TstVAVisManager;
    visManager->Initialize();
  #endif

  // UserAction classes
  runManager->SetUserAction(new TstVARunAction);
  runManager->SetUserAction(new TstVAPrimaryGeneratorAction);
  runManager->SetUserAction(new TstVASteppingAction);
  runManager->SetUserAction(new TstVAEventAction);

  // User interactions
  // Define (G)UI for interactive mode
 
  // User interactions
    G4UImanager * UI = G4UImanager::GetUIpointer();
  
  if(argc==1)
  { 
    // G4UIterminal is a (dumb) terminal.
    G4UIsession * session = new G4UIterminal;
    UI->ApplyCommand("/control/execute prerun.g4mac"); 
    session->SessionStart();
    delete session;
  }
  else
  // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }
  
  #ifdef G4VIS_USE
    delete visManager;
  #endif
  delete runManager;
  return 0;
}

