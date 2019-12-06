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
/// \file hadronic/Hadr00/Hadr00.cc
/// \brief Main program of the hadronic/Hadr00 example
//
//
//
// -------------------------------------------------------------
//      GEANT4 Hadr00
//
//  Application demonstrating Geant4 hadronic cross sections
//
//  Author: V.Ivanchenko 20 June 2008
//
//  Modified: 
//
// -------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "ActionInitialization.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) { ui = new G4UIExecutive(argc,argv); }

#ifdef G4MULTITHREADED  
  G4MTRunManager * runManager = new G4MTRunManager(); 

  // Number of threads can be defined via 3rd argument
  if (argc==4) {
    G4int nThreads = G4UIcommand::ConvertToInt(argv[3]);
    runManager->SetNumberOfThreads(nThreads);
  }
  G4cout << "##### Hadr00 started for " << runManager->GetNumberOfThreads()
         << " threads" << " #####" << G4endl;
#else
  G4RunManager * runManager = new G4RunManager();
  G4cout << "##### Hadr00 started in sequential mode"
         << " #####" << G4endl;
#endif

  //set mandatory initialization classes
  DetectorConstruction* det = new DetectorConstruction();
  runManager->SetUserInitialization(det);

  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = nullptr;
  G4String physName = "";

  //Physics List name defined via 3nd argument
  if (argc>=3) { physName = argv[2]; }

  //Physics List is defined via environment variable PHYSLIST
  if ("" == physName) {
    char* path = std::getenv("PHYSLIST");
    if (path) { physName = G4String(path); }
  }

  //if name is not known to the factory use FTFP_BERT
  if("" == physName || !factory.IsReferencePhysList(physName)) {
    physName = "FTFP_BERT";
  }

  //reference PhysicsList via its name
  phys = factory.GetReferencePhysList(physName);

  runManager->SetUserInitialization(phys);
  det->SetPhysicsList(phys);

  //set user action classes
  runManager->SetUserInitialization(new ActionInitialization(det));

  //initialize visualization
  G4VisManager* visManager = nullptr;

  //get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if (ui)  {
    //interactive mode
    visManager = new G4VisExecutive;
    visManager->Initialize();
    ui->SessionStart();
    delete ui;
  } else  {
    //batch mode  
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  //job termination 
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
