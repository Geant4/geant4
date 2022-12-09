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
/// \file Hadr08.cc
/// \brief Main program of the hadronic/Hadr08 example
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"
#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "FTFP_INCLXX.hh"
#include "G4UImanager.hh"
#include "G4RunManagerFactory.hh"
#include "ActionInitialization.hh"
#include "G4UIExecutive.hh"
#include "G4GenericBiasingPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main( int argc, char** argv ) {

  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  auto* runManager = G4RunManagerFactory::CreateRunManager();
  runManager->SetNumberOfThreads(4);

  G4VUserDetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization( detector );

  FTFP_BERT* physicsList = new FTFP_BERT;
  //FTFP_INCLXX* physicsList = new FTFP_INCLXX;

  G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics;
  biasingPhysics->Bias( "proton" );
  biasingPhysics->Bias( "neutron" );
  biasingPhysics->Bias( "pi+" );
  biasingPhysics->Bias( "pi-" );
  physicsList->RegisterPhysics( biasingPhysics );

  runManager->SetUserInitialization( physicsList );
  runManager->SetUserInitialization( new ActionInitialization );

  if ( ui ) {
    ui->SessionStart();
    delete ui;
  } else {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    G4UImanager * UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand(command+fileName);
  }

  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
