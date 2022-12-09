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
// (adapted from exampleB1)
// Author: A.Knaian (ara@nklabs.com), N.MacFadden (natemacfadden@gmail.com)


#include "G4MTRunManager.hh"
#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "G4ScoringManager.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

#include "QGSP_BIC.hh"
#include "G4StepLimiterPhysics.hh"

#include "FADetectorConstruction.hh"
#include "FAActionInitialization.hh"

int main(int argc,char** argv)
{
// Construct the default run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  G4int nThreads = 4;
  runManager->SetNumberOfThreads(nThreads);

	// Detect interactive mode (if no arguments) and define UI session
	//
	G4UIExecutive* ui = 0;
	if ( argc == 1 ) {
		ui = new G4UIExecutive(argc, argv);
	}

	// Set mandatory initialization classes
	//
	//G4ScoringManager* scoringManager =
   G4ScoringManager::GetScoringManager();

	// Detector construction initialization
	runManager->SetUserInitialization(new DetectorConstruction());

	// Physics list
	G4VModularPhysicsList* physicsList = new QGSP_BIC;
	physicsList->SetVerboseLevel(1);

	// Physics list - step limiter
	G4StepLimiterPhysics* stepLimiter = new G4StepLimiterPhysics;
	stepLimiter->SetApplyToAll(true);	// apply step limit to all particles. Default we set limit to DBL_MAX
	physicsList->RegisterPhysics(stepLimiter);

	// Physics list initialization
	runManager->SetUserInitialization(physicsList);
		
	// User action initialization
	runManager->SetUserInitialization(new ActionInitialization());
	
	// Initialize visualization
	//
	G4VisManager* visManager = new G4VisExecutive;
	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
	// G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();

	// Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();

	// Process macro or start UI session
	//
	if ( ! ui ) { 
		// batch mode
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command+fileName);
	}
	else { 
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		ui->SessionStart();
		delete ui;
	}

	// Job termination
	// Free the store: user actions, physics_list and detector_description are
	// owned and deleted by the run manager, so they should not be deleted 
	// in the main() program !
	
	delete visManager;
	delete runManager;
	//delete scoringManager;
}

