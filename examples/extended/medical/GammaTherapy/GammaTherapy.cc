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
//      GEANT4 ibrem
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- ibrem ----------------------------
//
//  Modified: 05.04.01 Vladimir Ivanchenko new design of ibrem
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#include "VisManager.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "RunAction.hh"
#include "Histo.hh"

#include <string>
#include <fstream>
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int main(int argc,char** argv) {

  G4int verbose = 1;
  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);

  //  G4String fName = argv[1];
  //  (Histo::GetInstance)->LoadParameters(fName);

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager();

  // set mandatory initialization classes
  DetectorConstruction* det = new DetectorConstruction();
  runManager->SetUserInitialization(det);

  runManager->SetUserInitialization(new PhysicsList());

  // visualization manager
  G4VisManager* visManager = new VisManager();
  visManager->Initialize();

  // set user action classes
  runManager->SetUserAction(new PrimaryGeneratorAction(det));
  runManager->SetUserAction(new RunAction());
  runManager->SetUserAction(new EventAction());
  runManager->SetUserAction(new TrackingAction());

  // get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(1 < verbose) UI->ListCommands("/testem/");

  if (argc==1)   // Define UI terminal for interactive mode
    {
     G4UIsession * session = new G4UIterminal;
#ifdef G4UI_USE_TCSH
     session = new G4UIterminal(new G4UItcsh);
#else
     session = new G4UIterminal();
#endif
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

  delete visManager;
  delete runManager;

  return 0;
}


