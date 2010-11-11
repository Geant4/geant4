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
// $Id: GammaRayTel.cc,v 1.19 2010-11-11 17:42:50 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelPhysicsList.hh"
#include "GammaRayTelPrimaryGeneratorAction.hh"
#include "GammaRayTelRunAction.hh"
#include "GammaRayTelEventAction.hh"

#include "QGSP_BIC.hh"

#ifdef G4ANALYSIS_USE
#include "GammaRayTelAnalysis.hh"
#endif


/* This global file is used to store relevant data for
   analysis with external tools */
std::ofstream outFile;

// This is the main function 
int main(int argc, char** argv)
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;
  
  // Set mandatory user initialization classes
  GammaRayTelDetectorConstruction* detector = 
    new GammaRayTelDetectorConstruction;
  runManager->SetUserInitialization(detector);

  // POSSIBILITY TO SELECT ANOTHER PHYSICS LIST
  //  do not use   GammaRayTelPhysicsList, this is old style and crashes at 
  //    program exit   
  //runManager->SetUserInitialization(new GammaRayTelPhysicsList);

  runManager->SetUserInitialization(new QGSP_BIC);

  // Set mandatory user action classes
  runManager->SetUserAction(new GammaRayTelPrimaryGeneratorAction);

#ifdef G4ANALYSIS_USE
  // Creation of the analysis manager
  GammaRayTelAnalysis* analysis = GammaRayTelAnalysis::getInstance();
#endif

  // Set optional user action classes
  GammaRayTelEventAction* eventAction = new GammaRayTelEventAction();
  GammaRayTelRunAction* runAction = new GammaRayTelRunAction();
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(runAction);

  // Set visualization and user interface
#ifdef G4VIS_USE
  // Visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  
  // Initialize G4 kernel
  //  runManager->Initialize();
  
  // Get the pointer to the UI manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  G4UIExecutive* ui = 0;
#ifdef G4UI_USE
  ui = new G4UIExecutive(argc, argv);
#endif

  if (ui)
    {
      /* prerunGammaRayTel.mac is loaded by default
	 unless a macro file is passed as the argument
	 of the executable */

      if(argc>1)
	{
	  G4String command = "/control/execute ";
	  for (int i=2; i<=argc; i++) 
	    {
	      G4String macroFileName = argv[i-1];
	      UImanager->ApplyCommand(command+macroFileName);
	    }
	}
      else  
	{
	  UImanager->ApplyCommand("/control/execute prerunGammaRayTel.mac");
	  ui->SessionStart();
	}
      delete ui;
    }

  // Job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
#ifdef G4ANALYSIS_USE
  delete analysis;
#endif
  delete runManager;
  return 0;
}
