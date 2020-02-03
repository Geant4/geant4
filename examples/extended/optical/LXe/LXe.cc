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
<<<<<<< HEAD
// $Id: LXe.cc 77782 2013-11-28 08:12:12Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
/// \file optical/LXe/LXe.cc
/// \brief Main program of the optical/LXe example
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4String.hh"

#include "LXePhysicsList.hh"
#include "LXeDetectorConstruction.hh"
#include "LXeActionInitialization.hh"

#include "LXeRecorderBase.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
<<<<<<< HEAD
=======
  //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) { ui = new G4UIExecutive(argc,argv); }

>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
#else
  G4RunManager * runManager = new G4RunManager;
#endif
<<<<<<< HEAD

  runManager->SetUserInitialization(new LXeDetectorConstruction());
  runManager->SetUserInitialization(new LXePhysicsList());
=======
  LXeDetectorConstruction* det = new LXeDetectorConstruction();
  runManager->SetUserInitialization(det);
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

  LXeRecorderBase* recorder = NULL; //No recording is done in this example

  runManager->SetUserInitialization(new LXeActionInitialization(recorder));

<<<<<<< HEAD
#ifdef G4VIS_USE
=======
  physicsList->RegisterPhysics(opticalPhysics);
  runManager->SetUserInitialization(physicsList);

  runManager->SetUserInitialization(new LXeActionInitialization(det));

  //initialize visualization
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // runManager->Initialize();
 
  // get the pointer to the UI manager and set verbosities
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

<<<<<<< HEAD
  if(argc==1){
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute vis.mac");
#endif
    if (ui->IsGUI())
       UImanager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();
    delete ui;
#endif
  }
  else{
=======
  if (ui) {
    //interactive mode
    UImanager->ApplyCommand("/control/execute vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  } else {
    //batch mode  
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
    G4String command = "/control/execute ";
    G4String filename = argv[1];
    UImanager->ApplyCommand(command+filename);
  }

//  if(recorder)delete recorder;

#ifdef G4VIS_USE
  delete visManager;
#endif

  // job termination
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
