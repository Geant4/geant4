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
//      GEANT4 hTest
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTest ----------------------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "hTestVisManager.hh"
#endif

#include "hTestDetectorConstruction.hh"
#include "hTestPhysicsList.hh"
#include "hTestPrimaryGeneratorAction.hh"
#include "hTestEventAction.hh"
#include "hTestTrackingAction.hh"
#include "hTestRunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int main(int argc,char** argv) {

  G4int verbose = 1;
  //choose the Random engine
  //  HepRandom::setTheEngine(new RanecuEngine);
  
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager();

  // set mandatory initialization classes
  hTestDetectorConstruction* det = new hTestDetectorConstruction();
  runManager->SetUserInitialization(det);
  det->SetVerbose(verbose);
  if(verbose >0) G4cout << "Detector Construction is defined" << G4endl;
  
  runManager->SetUserInitialization(new hTestPhysicsList(det));
  if(verbose >0) G4cout << "Physics List is defined" << G4endl;

#ifdef G4VIS_USE
  G4cout << "VisManager will be inicialized" << G4endl;
  // visualization manager
  G4VisManager* visManager = new hTestVisManager;
  visManager->Initialize();
#endif 

  G4cout << "User actions will be initialized" << G4endl;
 
  // set user action classes
  runManager->SetUserAction(new hTestPrimaryGeneratorAction(det));
  if(verbose >0) G4cout << "hTestPrimaryGeneratorAction is defined" << G4endl;

  runManager->SetUserAction(new hTestRunAction());
  if(verbose >0) G4cout << "hTestRunAction is defined" << G4endl;

  hTestEventAction* event = new hTestEventAction(det);
  runManager->SetUserAction(event);
  if(verbose >0) G4cout << "EventAction is defined" << G4endl;

  runManager->SetUserAction(new hTestTrackingAction());
  if(verbose >0) G4cout << "TrackingAction is defined" << G4endl;
  
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  if(1 < verbose) UI->ListCommands("/hTest/");

  if (argc==1)   // Define UI terminal for interactive mode  
    { 
     // Initialize G4 kernel
     runManager->Initialize();

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

     // Initialize G4 kernel
     G4int nev = det->GetNumberOfEvents();
     if(nev > 0) {
       G4cout << "Start initialisation for hTest" << G4endl;
       runManager->Initialize();

       if(verbose >0) {
         G4cout << "Start event loop for " << nev  
                << " events" << G4endl;
       }

       runManager->BeamOn(nev);
     }

     // next file
     if(argc==3) {
       if(verbose >0) G4cout << "Second mac file is applied" << G4endl; 
       UI->ApplyCommand(command+argv[2]);
     }
    }
    
  // job termination
  
#ifdef G4VIS_USE
  delete visManager;
#endif

  G4cout << "runManager will be deleted" << G4endl;  
  delete runManager;
  G4cout << "runManager is deleted" << G4endl;  

  return 0;
}

