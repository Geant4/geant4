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
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo
//
//   *******************************************************
//   *                  Ultra.cc
//   *******************************************************


#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "UltraRunAction.hh"
#include "UltraDetectorConstruction.hh"
#include "UltraPrimaryGeneratorAction.hh"
#include "UltraPhysicsList.hh"
#include "UltraEventAction.hh"
//#include "CLHEP/Random/RanluxEngine.h"

int main(int argc,char** argv) {

//choose the Random engine from CLHEP 
//(lets use C++ implementation of Jame's RANLUX generator)
 
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanluxEngine);
  G4RunManager* runManager = new G4RunManager;

  // UserInitialization classes - mandatory
  UltraDetectorConstruction* detector = new UltraDetectorConstruction;
  UltraPhysicsList* list = new UltraPhysicsList();
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(list);

  // UserAction classes - optional
  UltraPrimaryGeneratorAction* PrimGenAct = new UltraPrimaryGeneratorAction();
  runManager->SetUserAction(PrimGenAct);
  UltraRunAction* RunAct = new UltraRunAction();
  runManager->SetUserAction(RunAct);
  UltraEventAction* EvAct = new UltraEventAction(RunAct);
  runManager->SetUserAction(EvAct);

#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  //Initialize G4 kernel
  runManager->Initialize();

  // Get the Pointer to the UI Manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // User interactions
  // Define (G)UI for interactive mode
  if(argc==1)
  {
    // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_XM
    G4UIsession* session = new G4UIXm(argc,argv);
#else
    G4UIsession* session = new G4UIterminal();
#endif
    UImanager->ApplyCommand("/control/execute Visualisation.mac");
    session->SessionStart();
    delete session;
  }
  else
  // Batch mode
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete  visManager;
#endif

  delete runManager;

  return 0;
}

