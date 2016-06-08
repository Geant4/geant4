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
// $Id: PersistentEx01.cc,v 1.5.4.2 2001/06/28 20:18:46 gunter Exp $
// GEANT4 tag $Name:  $
//
//
// --------------------------------------------------------------
//      GEANT4 - PersistentEx01
//
// --------------------------------------------------------------
// Comments
//
//
// --------------------------------------------------------------


#include "PersEx01DetectorConstruction.hh"
#include "PersEx01PhysicsList.hh"
#include "PersEx01RunAction.hh"
#include "PersEx01PrimaryGeneratorAction.hh"
#include "PersEx01EventAction.hh"
#include "PersEx01StackingAction.hh"

#include "G4ios.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4PersistencyManager.hh"

int main(int argc, char** argv)
{
  // Persistency Manager
  G4PersistencyManager * persistencyManager =
                          G4PersistencyManager::GetPersistencyManager();

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // Detector geometry
  runManager->SetUserInitialization(new PersEx01DetectorConstruction);
  runManager->SetUserInitialization(new PersEx01PhysicsList);
  
  // UserAction classes.
  runManager->SetUserAction(new PersEx01RunAction);
  runManager->SetUserAction(new PersEx01PrimaryGeneratorAction);
  runManager->SetUserAction(new PersEx01EventAction);
  runManager->SetUserAction(new PersEx01StackingAction);

  // get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if(argc==1)
  // Define (G)UI terminal for interactive mode
  {
    G4UIsession* session = new G4UIterminal;
    UI->ApplyCommand("/control/execute prerunPEx01.mac");
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

  // job termination
  delete runManager;
  delete persistencyManager;

  return 0;
}

