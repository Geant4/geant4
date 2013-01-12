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
// $Id: $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "StackingAction.hh"

#ifdef G4_USE_ROOT
#include "ROOTAnalysis.hh"
#endif

#include "G4Timer.hh"
#include "G4VisExecutive.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  G4Timer* theTimer = new G4Timer();

  theTimer->Start();

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);
  PhysicsList* physicsList = new PhysicsList();
  runManager->SetUserInitialization(physicsList);
  
  // primary generator
  PrimaryGeneratorAction* primary = new PrimaryGeneratorAction(detector);
  runManager->SetUserAction(primary);
        
  // set user action classes
  RunAction*      runAct = new RunAction(physicsList);
  EventAction*    evtAct = new EventAction(runAct);
  SteppingAction* stpAct = new SteppingAction();
  
  runManager->SetUserAction(runAct);
  runManager->SetUserAction(evtAct);
  runManager->SetUserAction(stpAct);
  runManager->SetUserAction(new StackingAction());

  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc!=1)   // batch mode  
    {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
    }
    
  else           //define visualization and UI terminal for interactive mode
    { 
   G4VisManager* visManager = new G4VisExecutive;
   visManager->Initialize();   
     
#ifdef _WIN32
    G4UIsession * session = new G4UIterminal();
#else
    G4UIsession * session = new G4UIterminal(new G4UItcsh);
#endif
     session->SessionStart();
     delete session;
     
     delete visManager;
    }

  //Close manually the ROOT file here, if it is the case
#ifdef G4_USE_ROOT
  G4int returnCode = ROOTAnalysis::getInstance()->CheckPhysicalResults();
  if (!returnCode)
    G4cout << "TEST PASSED" << G4endl;
  G4cout << "Going to close the ROOT file....";
  ROOTAnalysis::getInstance()->CloseFile();
  G4cout << "done" << G4endl;

#endif

  // job termination
  //   
  delete runManager;
  theTimer->Stop();

  G4cout << "Execution terminated" << G4endl;
  G4cout << (*theTimer) << G4endl;
  delete theTimer;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
