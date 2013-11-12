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

#include "Tst69DetectorConstruction.hh"
#include "Tst69ActionInitialization.hh"
#include "Tst69MacroLocation.hh"

#include "G4UImanager.hh"
#include "G4UIExecutive.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif
#include "G4ios.hh"

#include "G4PhysListFactory.hh"

#include <cstdlib>

int main(int argc,char** argv) {

  // Set the default random engine to RanecuEngine
  CLHEP::RanecuEngine defaultEngine;
  G4Random::setTheEngine(&defaultEngine);

#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Run manager
#ifdef G4MULTITHREADED
  char env[]="G4FORCENUMBEROFTHREADS=4";
  putenv(env);

  G4MTRunManager* runManager = new G4MTRunManager;
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // UserInitialization classes
  runManager->SetUserInitialization(new Tst69DetectorConstruction);

  // Select the physics list using the TEST69_PHYSLIST environment variable
  const char *physListEnvVar = getenv("TEST69_PHYSLIST");
  G4String physicsListName;
  if(physListEnvVar == NULL)
    physicsListName = "QGSP_INCLXX";
  else
    physicsListName = physListEnvVar;
  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = factory.GetReferencePhysList(physicsListName);
  runManager->SetUserInitialization(phys);

  // UserAction classes
  runManager->SetUserInitialization(new Tst69ActionInitialization);

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if(argc==1)
  {
    // Define the UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macroLocation+"/vis.mac"); 
#endif
    ui->SessionStart();
    delete ui;
#endif
  }
  else
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

