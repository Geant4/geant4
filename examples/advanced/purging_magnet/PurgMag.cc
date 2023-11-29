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
// Code developed by:
//  S.Larsson
//
//    *********************
//    *                   *
//    *    PurgMag.cc     *
//    *                   *
//    *********************
//
//
// Comments: Main program for the Purgin Magnet example.
//

#include "G4Types.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4AnalysisManager.hh"
#include "PurgMagDetectorConstruction.hh"
#include "PurgMagPhysicsList.hh"
#include "PurgMagActionInitializer.hh"

int main(int argc,char** argv) 
{
 // Construct the default run manager
 //
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  G4int nThreads = 4;
  runManager->SetNumberOfThreads(nThreads);

  // set mandatory initialization classes
  runManager->SetUserInitialization(new PurgMagDetectorConstruction);
  runManager->SetUserInitialization(new PurgMagPhysicsList);
  runManager->SetUserInitialization(new PurgMagActionInitializer());

  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  //Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();


  if (argc==1)   // Define UI session for interactive mode.
    {
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      UImanager->ApplyCommand("/control/execute vis.mac");
      ui->SessionStart();
      delete ui;
    }
  else           // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

  // Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();

  // job termination
  delete visManager;

  delete runManager;

  return 0;
}

