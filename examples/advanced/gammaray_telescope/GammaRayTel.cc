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
//
// ------------------------------------------------------------
//      GEANT 4 main program
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTel example main program ------
//           by F.Longo, R.Giannitrapani & G.Santin (29 nov 2000)
//           See README file for details on this example
//  20.11.01 G.Santin: new analysis management, and some modification in the
//                     construction of some Action's
// ************************************************************

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#include "G4VisExecutive.hh"

#include "G4UIExecutive.hh"

#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelPhysicsList.hh"
#include "GammaRayTelActionInitializer.hh"

//#include "QGSP_BIC.hh"
#include "FTFP_BERT.hh"

#include "GammaRayTelAnalysis.hh"

// This is the main function
int main(int argc, char** argv)
{
  // Construct the default run manager
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  //  runManager->SetNumberOfThreads(1);
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // Set mandatory user initialization classes
  GammaRayTelDetectorConstruction* detector =
    new GammaRayTelDetectorConstruction;
  runManager->SetUserInitialization(detector);

  // POSSIBILITY TO SELECT ANOTHER PHYSICS LIST
  //  do not use   GammaRayTelPhysicsList, this is old style and crashes at
  //    program exit
  runManager->SetUserInitialization(new GammaRayTelPhysicsList);

  //  runManager->SetUserInitialization(new QGSP_BIC);
  //runManager->SetUserInitialization(new FTFP_BERT);

  //Initialize actions
  runManager->SetUserInitialization(new GammaRayTelActionInitializer());

  // Creation of the analysis manager
  GammaRayTelAnalysis* analysis = GammaRayTelAnalysis::getInstance();

  // Set visualization and user interface
  // Visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // Initialize G4 kernel
  //  runManager->Initialize();

  // Get the pointer to the UI manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  else
    {
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      if (ui->IsGUI())
	{
	  /* prerunGammaRayTel.mac is loaded by default */
	  UImanager->ApplyCommand("/control/execute prerunGammaRayTel.mac");
	  ui->SessionStart();
	}
      delete ui;
    }
  // Job termination
  delete visManager;
  delete analysis;
  delete runManager;
  return 0;
}
