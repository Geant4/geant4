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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"

#include "SAXSDetectorConstruction.hh"
#include "SAXSPhysicsList.hh"
#include "SAXSActionInitialization.hh"
#include "SAXSRunAction.hh"
#include "SAXSEventAction.hh"
#include "SAXSSteppingAction.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4Timer.hh"

#include "Randomize.hh"
#include <time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{

  G4Timer* theTimer = new G4Timer();
  theTimer->Start();

  //Choose the Random engine and reset the seed (before the RunManager)
  /*
    G4int seed = time(0);
    G4cout << "Random seed: " << seed << G4endl;
    CLHEP::HepRandom::setTheSeed(seed);
  */
  //Construct the run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  int vNumberOfThreads = 2;
  if (argc>2) {
    vNumberOfThreads = atoi(argv[2]);
  }
  if (vNumberOfThreads > 0) {
    runManager->SetNumberOfThreads(vNumberOfThreads);
  }

  //Set mandatory initialization classes
  runManager->SetUserInitialization(new SAXSDetectorConstruction);
  runManager->SetUserInitialization(new SAXSPhysicsList());

  //Set user action classes
  runManager->SetUserInitialization(new SAXSActionInitialization());
 
  //Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  //No arguments: set the batch mode
  if (argc!=1) {
    //Batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  } else {
    //Visualization manager
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    //Define UI session for interactive mode
    G4UIExecutive* ui = new G4UIExecutive(argc,argv);
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI())
      UImanager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();

    delete ui;
    delete visManager;
  }

  //Job termination
  delete runManager;

  theTimer->Stop();
  G4cout << "Execution completed" << G4endl;
  G4cout << (*theTimer) << G4endl;
  delete theTimer;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

