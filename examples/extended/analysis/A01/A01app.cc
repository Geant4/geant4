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
/// \file analysis/A01/A01app.cc
/// \brief Main program of the analysis/A01 example
//
// $Id$
// --------------------------------------------------------------
//
// --------------------------------------------------------------
//      GEANT 4 - A01app
// --------------------------------------------------------------
// Comments
//   Tutorial for Geant4 lectures
// --------------------------------------------------------------


#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "A01DetectorConstruction.hh"
#include "A01PhysicsList.hh"
#include "A01PrimaryGeneratorAction.hh"

#include "A01EventAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv)
{
  // RunManager construction
  G4RunManager* runManager = new G4RunManager;

#ifdef G4VIS_USE
  // Visualization manager construction
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // mandatory user initialization classes
  runManager->SetUserInitialization(new A01DetectorConstruction);
  runManager->SetUserInitialization(new A01PhysicsList);

  // initialize Geant4 kernel
  runManager->Initialize();

  // mandatory user action class
  runManager->SetUserAction(new A01PrimaryGeneratorAction);

  // optional user action classes
  runManager->SetUserAction(new A01EventAction);

  if(argc>1)
  // execute an argument macro file if exist
  {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else
  // start interactive session
  {
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    ui->SessionStart();
    delete ui;
#endif
  }

#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;

  return 0;
}

