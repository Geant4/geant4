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
// -------------------------------------------------------------
//      GEANT4 test31
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- test31 ----------------------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of test31 
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "VisManager.hh"
#endif

#include "test31DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "test31PrimaryGeneratorAction.hh"
#include "test31EventAction.hh"
#include "test31TrackingAction.hh"
#include "test31RunAction.hh"

#include "G4Timer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int main(int argc,char** argv) {

  G4Timer* timer = new G4Timer();
  timer->Start();

  G4int verbose = 1;
  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager();

  // set mandatory initialization classes
  test31DetectorConstruction* det = new test31DetectorConstruction();
  runManager->SetUserInitialization(det);
  det->SetVerbose(verbose);
  if(verbose >0) G4cout << "Detector Construction is defined" << G4endl;
  
  runManager->SetUserInitialization(new PhysicsList());
  if(verbose >0) G4cout << "Physics List is defined" << G4endl;

#ifdef G4VIS_USE
  G4cout << "VisManager will be inicialized" << G4endl;
  // visualization manager
  G4VisManager* visManager = new VisManager;
  visManager->Initialize();
#endif 

  G4cout << "User actions will be initialized" << G4endl;
 
  // set user action classes
  runManager->SetUserAction(new test31PrimaryGeneratorAction(det));
  if(verbose >0) G4cout << "test31PrimaryGeneratorAction is defined" << G4endl;

  runManager->SetUserAction(new test31RunAction());
  if(verbose >0) G4cout << "test31RunAction is defined" << G4endl;

  test31EventAction* event = new test31EventAction(det);
  runManager->SetUserAction(event);
  if(verbose >0) G4cout << "EventAction is defined" << G4endl;

  runManager->SetUserAction(new test31TrackingAction());
  if(verbose >0) G4cout << "TrackingAction is defined" << G4endl;

  // get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(1 < verbose) UI->ListCommands("/test31/");
     
  if (argc==1)   // Define UI terminal for interactive mode
    {

     G4UIsession * session = new G4UIterminal;
     UI->ApplyCommand("/control/execute init.mac");
     session->SessionStart();
     delete session;
    }
  else if (argc>1) // Batch mode with 1 or more files
    {
     if(verbose >0) G4cout << "UI interface is started" << G4endl;
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
    }
    
 // job termination
  
#ifdef G4VIS_USE
  delete visManager;
#endif

  timer->Stop();
  G4cout << "  "  << *timer << G4endl;
  delete timer;

  //  G4cout << "runManager will be deleted" << G4endl;  
  delete runManager;
  return 0;
}

