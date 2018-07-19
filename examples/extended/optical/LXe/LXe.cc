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
// $Id: LXe.cc 110280 2018-05-17 14:50:16Z gcosmo $
//
/// \file optical/LXe/LXe.cc
/// \brief Main program of the optical/LXe example
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4String.hh"

#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"

#include "LXeDetectorConstruction.hh"

#include "LXeActionInitialization.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{

  //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc,argv);

#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  G4int nThreads = std::min(G4Threading::G4GetNumberOfCores(), 4);
  runManager->SetNumberOfThreads(nThreads);
  G4cout << "===== LXe is started with "
         <<  runManager->GetNumberOfThreads() << " threads =====" << G4endl;
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  runManager->SetUserInitialization(new LXeDetectorConstruction());

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  opticalPhysics->SetWLSTimeProfile("delta");

  opticalPhysics->SetScintillationYieldFactor(1.0);
  opticalPhysics->SetScintillationExcitationRatio(0.0);

  opticalPhysics->SetMaxNumPhotonsPerStep(100);
  opticalPhysics->SetMaxBetaChangePerStep(10.0);

  opticalPhysics->SetTrackSecondariesFirst(kCerenkov, true);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation, true);

  physicsList->RegisterPhysics(opticalPhysics);
  runManager->SetUserInitialization(physicsList);

  runManager->SetUserInitialization(new LXeActionInitialization());

  //initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  //get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (ui)  {
    //interactive mode
    UImanager->ApplyCommand("/control/execute vis.mac");
    if (ui->IsGUI()) {
     UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }
  else  {
    //batch mode  
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  // job termination
  delete visManager;
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
