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
// $Id: exampleN03.cc,v 1.36 2008/11/24 13:54:38 lgarnier Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "ExN03DetectorConstruction.hh"
#include "ExN03PhysicsList.hh"
#include "ExN03PrimaryGeneratorAction.hh"
#include "ExN03RunAction.hh"
#include "ExN03EventAction.hh"
#include "ExN03SteppingAction.hh"
#include "ExN03SteppingVerbose.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#if defined(G4UI_USE_TCSH)
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#elif defined(G4UI_USE_XM)
#include "G4UIXm.hh"
#elif defined(G4UI_USE_WIN32)
#include "G4UIWin32.hh"
#elif defined(G4UI_USE_QT)
#include "G4UIQt.hh"
#include "G4Qt.hh"
#else
#include "G4UIterminal.hh"
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // User Verbose output class
  //
  G4VSteppingVerbose::SetInstance(new ExN03SteppingVerbose);
     
  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  //
  ExN03DetectorConstruction* detector = new ExN03DetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new ExN03PhysicsList;
  runManager->SetUserInitialization(physics);
    
  // Set user action classes
  //
  G4VUserPrimaryGeneratorAction* gen_action = 
                          new ExN03PrimaryGeneratorAction(detector);
  runManager->SetUserAction(gen_action);
  //
  ExN03RunAction* run_action = new ExN03RunAction;  
  runManager->SetUserAction(run_action);
  //
  ExN03EventAction* event_action = new ExN03EventAction(run_action);
  runManager->SetUserAction(event_action);
  //
  G4UserSteppingAction* stepping_action =
                    new ExN03SteppingAction(detector, event_action);
  runManager->SetUserAction(stepping_action);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
  
#ifdef G4VIS_USE
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UI = G4UImanager::GetUIpointer();      
  
  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);    
    }
  else           // interactive mode : define visualization UI terminal
    {
      G4UIsession* session = 0;
#if defined(G4UI_USE_TCSH)
      session = new G4UIterminal(new G4UItcsh);      
#elif defined(G4UI_USE_XM)
      session = new G4UIXm(argc,argv);
      UI->ApplyCommand("/control/execute visTutor/gui.mac");      
#elif defined(G4UI_USE_WIN32)
      session = new G4UIWin32();
      UI->ApplyCommand("/control/execute visTutor/gui.mac");      
#elif defined(G4UI_USE_QT)
      session = new G4UIQt(argc,argv);
      UI->ApplyCommand("/control/execute visTutor/gui.mac");      
#else
      session = new G4UIterminal();
#endif
#ifdef G4VIS_USE
      UI->ApplyCommand("/control/execute vis.mac");     
#endif
      session->SessionStart();
      delete session;
    }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
#ifdef G4VIS_USE
  delete visManager;
#endif                
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
