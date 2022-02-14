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
<<<<<<< HEAD
// $Id: Hadr00.cc 81073 2014-05-20 10:23:13Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
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

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "ActionInitialization.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  //choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine());

  auto* runManager = G4RunManagerFactory::CreateRunManager();

  // Number of threads can be defined via 3rd argument
  if (argc==4) {
    G4int nThreads = G4UIcommand::ConvertToInt(argv[3]);
    runManager->SetNumberOfThreads(nThreads);
  }
  G4cout << "##### Hadr00 started for " << runManager->GetNumberOfThreads() 
         << " threads" << " #####" << G4endl;

  //set mandatory initialization classes
  DetectorConstruction* det = new DetectorConstruction();
  runManager->SetUserInitialization(det);

  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = 0;
  G4String physName = "";

  // Physics List name defined via 3nd argument
  if (argc>=3) { physName = argv[2]; }

<<<<<<< HEAD
  // Physics List is defined via environment variable PHYSLIST
  if("" == physName) {
    char* path = getenv("PHYSLIST");
=======
  //Physics List is defined via environment variable PHYSLIST
  if ("" == physName) {
    char* path = std::getenv("PHYSLIST");
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
    if (path) { physName = G4String(path); }
  }

  // if name is not known to the factory use FTFP_BERT
  if("" == physName || !factory.IsReferencePhysList(physName)) {
    physName = "FTFP_BERT"; 
  }

  // reference PhysicsList via its name
  phys = factory.GetReferencePhysList(physName);

  runManager->SetUserInitialization(phys);

  //set user action classes
  runManager->SetUserInitialization(new ActionInitialization(det));

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
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
