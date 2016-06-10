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
// $Id: $
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


#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4UImanager.hh"

#include "RE01DetectorConstruction.hh"
#include "RE01PhysicsList.hh"
#include "QGSP_BERT.hh"

#include "RE01ActionInitialization.hh"


#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv)
{
#ifdef G4MULTITHREADED
 G4MTRunManager * runManager = new G4MTRunManager;
 //runManager->SetNumberOfThreads(4);
#else
 G4RunManager * runManager = new G4RunManager;
#endif

#ifdef G4VIS_USE
  // Visualization manager construction
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  runManager->SetUserInitialization(new RE01DetectorConstruction);
  runManager->SetUserInitialization(new RE01PhysicsList);
  //runManager->SetUserInitialization(new QGSP_BERT);
  
  runManager->SetUserInitialization(new RE01ActionInitialization);
  
  runManager->Initialize();
  
  if(argc==1)
  {
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    ui->SessionStart();
    delete ui;
#endif
  }
  else
  {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();  
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

