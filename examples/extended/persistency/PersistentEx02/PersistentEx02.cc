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
// $Id: PersistentEx02.cc,v 1.4 2001/07/11 09:58:14 gunter Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
//
// --------------------------------------------------------------
//      GEANT4 - PersistentEx02
//
// --------------------------------------------------------------
// Comments
//
//
// --------------------------------------------------------------


#include "PersEx02DetectorConstruction.hh"
#include "PersEx02PhysicsList.hh"
#include "PersEx02RunAction.hh"
#include "PersEx02PrimaryGeneratorAction.hh"
#include "PersEx02EventAction.hh"
#include "PersEx02StackingAction.hh"

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
  runManager->SetUserInitialization(new PersEx02DetectorConstruction);
  runManager->SetUserInitialization(new PersEx02PhysicsList);
  
  // UserAction classes.
  runManager->SetUserAction(new PersEx02RunAction);
  runManager->SetUserAction(new PersEx02PrimaryGeneratorAction);
  runManager->SetUserAction(new PersEx02EventAction);
  runManager->SetUserAction(new PersEx02StackingAction);

  // get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if(argc==1)
  // Define (G)UI terminal for interactive mode
  {
    G4UIsession* session = new G4UIterminal;
    UI->ApplyCommand("/control/execute prerunPEx02.mac");
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

