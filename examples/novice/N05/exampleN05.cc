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
// $Id: exampleN05.cc,v 1.10 2002-01-09 17:24:16 ranjard Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN05
//
// --------------------------------------------------------------
// Comments
//
// * Example of a main program making use of parameterisation
//   ie "Fast Simulation".
//-------------------------------------------------------------------

//------------------------------
// Run, Generator etc.. actions:
//------------------------------
#include "ExN05RunAction.hh"
#include "ExN05PrimaryGeneratorAction.hh"
#include "ExN05EventAction.hh"
#include "ExN05SteppingAction.hh"

//---------------------------
// Parameterisation manager:
//---------------------------
#include "G4GlobalFastSimulationManager.hh"

//--------------------
// Detector:
//--------------------
#include "ExN05DetectorConstruction.hh"

//----------------------------------
// ExN05PhysicsList makes use of the
// G4ParameterisationManagerProcess
//----------------------------------
#include "ExN05PhysicsList.hh"

#include "G4UIterminal.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"

#ifdef G4VIS_USE
#include "ExN05VisManager.hh"
#endif

#include "G4ios.hh"


int main(int argc, char** argv) {

  //-------------------------------
  // Initialization of Run manager
  //-------------------------------
  G4cout << "RunManager construction starting...." << G4endl;
  G4RunManager * runManager = new G4RunManager;

  // Detector geometry
  ExN05DetectorConstruction* Detector = new ExN05DetectorConstruction();
  runManager->SetUserInitialization(Detector);

  // PhysicsList: (including G4FastSimulationManagerProcess)
  runManager->SetUserInitialization(new ExN05PhysicsList);

  // UserAction classes.
  runManager->SetUserAction(new ExN05RunAction);
  runManager->SetUserAction(new ExN05PrimaryGeneratorAction);
  runManager->SetUserAction(new ExN05EventAction);
  runManager->SetUserAction(new ExN05SteppingAction);

  // Inizialize Run manager
  runManager->Initialize();

  // Close the "fast simulation": will
  // trigger the ghost geomtries construction:
  G4GlobalFastSimulationManager::GetGlobalFastSimulationManager()->
    CloseFastSimulation();

    
  //----------------
  // Visualization:
  //----------------
#ifdef G4VIS_USE
  G4cout << "Instantiating Visualization Manager......." << G4endl;
  G4VisManager* visManager = new ExN05VisManager;
  visManager -> Initialize ();
#endif

  //Setup commandes:
  G4UImanager * UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/Step/Verbose 0");
  UI->ApplyCommand("/tracking/Verbose 1");
  UI->ApplyCommand("/gun/particle e-");
  UI->ApplyCommand("/gun/energy 100 MeV");
  UI->ApplyCommand("/gun/direction 0 0 1");
  UI->ApplyCommand("/gun/position 0 0 0");
  UI->ApplyCommand("/gun/direction 0 .3 1.");

  if(argc==1)
  {
    //--------------------------
    // Define (G)UI
    //--------------------------
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

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  return EXIT_SUCCESS;
}
