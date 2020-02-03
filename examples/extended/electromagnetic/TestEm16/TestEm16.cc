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
/// \file electromagnetic/TestEm16/TestEm16.cc
/// \brief Main program of the electromagnetic/TestEm16 example
//
<<<<<<< HEAD
// $Id: TestEm16.cc 85257 2014-10-27 08:50:30Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "SteppingVerbose.hh"

#ifdef G4VIS_USE
 #include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

<<<<<<< HEAD
  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
=======
  //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc,argv);

  //construct the default run manager
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    G4int nThreads = G4Threading::G4GetNumberOfCores();
    if (argc==3) nThreads = G4UIcommand::ConvertToInt(argv[2]);
    runManager->SetNumberOfThreads(nThreads);
#else
    G4VSteppingVerbose::SetInstance(new SteppingVerbose);
    G4RunManager* runManager = new G4RunManager;
#endif

  //set mandatory initialization classes
  DetectorConstruction* det;
  runManager->SetUserInitialization(det = new DetectorConstruction);
  runManager->SetUserInitialization(new PhysicsList);

  //set user action classes
  runManager->SetUserInitialization(new ActionInitialization(det));

  //get the pointer to the User Interface manager
<<<<<<< HEAD
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (argc!=1)   // batch mode  
    {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
    }
    
  else           //define visualization and UI terminal for interactive mode
    {
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
=======
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (ui)  {
    //interactive mode
    visManager = new G4VisExecutive;
    visManager->Initialize();
    // define icons before SessionStart
    if (ui->IsGUI()) { UImanager->ApplyCommand("/control/execute gui.mac"); }
    ui->SessionStart();
    delete ui;
  } else  {
    //batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

#ifdef G4UI_USE
  G4UIExecutive * ui = new G4UIExecutive(argc,argv);      
  ui->SessionStart();
  delete ui;
#endif
          
#ifdef G4VIS_USE
  delete visManager;
#endif
    }

  //job termination
  //
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
