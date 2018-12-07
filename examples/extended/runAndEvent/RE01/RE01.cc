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
/// \file runAndEvent/RE01/RE01.cc
/// \brief Main program of the runAndEvent/RE01 example
//
//
//
//
// --------------------------------------------------------------
//      GEANT4 - RE01 exsample code
//
// --------------------------------------------------------------
// Comments
//
//
// --------------------------------------------------------------

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4UImanager.hh"

#include "RE01DetectorConstruction.hh"
#include "RE01CalorimeterROGeometry.hh"
#include "QGSP_BERT.hh"
#include "G4UnknownDecayPhysics.hh"
#include "G4ParallelWorldPhysics.hh"
#include "RE01ActionInitialization.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main(int argc,char** argv)
{
 // Instantiate G4UIExecutive if there are no arguments (interactive mode)
 G4UIExecutive* ui = nullptr;
 if ( argc == 1 ) {
   ui = new G4UIExecutive(argc, argv);
 }

#ifdef G4MULTITHREADED
 G4MTRunManager * runManager = new G4MTRunManager;
 //runManager->SetNumberOfThreads(4);
#else
 G4RunManager * runManager = new G4RunManager;
#endif

  // Visualization manager construction
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4String parallelWorldName = "ReadoutWorld";
  G4VUserDetectorConstruction* detector
     = new RE01DetectorConstruction();
  detector->RegisterParallelWorld(
       new RE01CalorimeterROGeometry(parallelWorldName));
  runManager->SetUserInitialization(detector);

  G4VModularPhysicsList* physicsList = new QGSP_BERT;
  physicsList->RegisterPhysics(new G4UnknownDecayPhysics);
  physicsList->RegisterPhysics(
       new G4ParallelWorldPhysics(parallelWorldName));
  runManager->SetUserInitialization(physicsList);

  runManager->SetUserInitialization(
       new RE01ActionInitialization);

  runManager->Initialize();

  if(ui)
  {
    ui->SessionStart();
    delete ui;
  }
  else
  {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  delete visManager;

  delete runManager;

  return 0;
}

