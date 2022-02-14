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
/// \file persistency/P01/exampleP01.cc
/// \brief Main program of the persistency/P01 example
//
//
<<<<<<< HEAD
// $Id: exampleP01.cc 82130 2014-06-11 09:26:44Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExP01DetectorConstruction.hh"
#include "ExP01PrimaryGeneratorAction.hh"
#include "ExP01RunAction.hh"
#include "ExP01EventAction.hh"
#include "ExP01SteppingAction.hh"
#include "ExP01SteppingVerbose.hh"

#include "FTFP_BERT.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new ExP01SteppingVerbose);

  // Run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);

  // UserInitialization classes (mandatory)
  ExP01DetectorConstruction* ExP01detector = new ExP01DetectorConstruction;
  runManager->SetUserInitialization(ExP01detector);

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  runManager->SetUserInitialization(physicsList);

  // Visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
   
  // UserAction classes
  runManager->SetUserAction(new ExP01PrimaryGeneratorAction(ExP01detector));
  runManager->SetUserAction(new ExP01RunAction);
  runManager->SetUserAction(new ExP01EventAction);
  runManager->SetUserAction(new ExP01SteppingAction);

  //Initialize G4 kernel
  runManager->Initialize();

  //get the pointer to the User Interface manager
  G4UImanager * UImanager = G4UImanager::GetUIpointer();

  if(ui)
  // Define (G)UI terminal for interactive mode
  {
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
#endif
  }
  else
  // Batch mode
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

