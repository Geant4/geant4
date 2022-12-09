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
// Author: Haegin Han
// Contributor: Min Cheol Han, Bangho Shin, Chansoo Choi, Yeon Soo Yeom, 
//              Jonghwi Jeong, Chan Hyeong Kim
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Geant4 Contributors: J. Allison and S. Guatelli
//
#include "TETDetectorConstruction.hh"
#include "TETModelImport.hh"
#include "TETActionInitialization.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "QGSP_BIC_HP.hh"

void PrintUsage(){
	G4cerr<< "Usage: ./External -m [MACRO] -o [OUTPUT] -f (option for MRCP-AF phantom)"  <<G4endl;
	G4cerr<< "Example: ./External -m run.mac -o run.out (-f)" <<G4endl;
}

int main(int argc,char** argv) 
{
	// Read the arguments for batch mode
	//
	G4String macro;
	G4String output;
	G4bool   isAF(false);
	G4UIExecutive* ui = nullptr;

	for ( G4int i=1; i<argc; i++ ) {
		// macro file name
		if ( G4String(argv[i]) == "-m" ) {
			macro = argv[i+1];
			i++;
		}
		// output file name
		else if ( G4String(argv[i]) == "-o" ) {
			output = argv[i+1];
			i++;
		}
		// switch for MRCP-AF phantom
		else if ( G4String(argv[i]) == "-f" ) {
			isAF = true;
		}
		else {
			PrintUsage();
			return 1;
		}
	}

	// print usage when there are more than six arguments
	if ( argc>6 ){
		PrintUsage();
		return 1;
	}

	// Detect interactive mode (if no macro file name) and define UI session
	//
	if ( !macro.size() ) {
		ui = new G4UIExecutive(argc, argv);
	}
	// default output file name
	else if ( !output.size() ) output = macro + ".out";

	// Choose the Random engine
	//
	//G4Random::setTheSeed(time(0));

	// Construct the default run manager
	//
       auto* runManager = G4RunManagerFactory::CreateRunManager();
       G4int nThreads = 4;
       runManager->SetNumberOfThreads(nThreads);
 
	// Set a class to import phantom data
	//
	auto* tetData = new TETModelImport(isAF, ui);

	// Set mandatory initialisation classes
	//
	// detector construction
	runManager->SetUserInitialization(new TETDetectorConstruction(tetData));
	// physics list
	//runManager->SetUserInitialization(new TETPhysicsList());
	runManager->SetUserInitialization(new QGSP_BIC_HP());
	// user action initialisation
	runManager->SetUserInitialization(new TETActionInitialization(tetData, output));
    
	// Visualization manager
	//
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialise();

	// Process macro or start UI session
	//
	G4UImanager* UImanager = G4UImanager::GetUIpointer();

	if ( ! ui ){
		// batch mode
		G4String command = "/control/execute ";
		UImanager->ApplyCommand(command+macro);
	}
	else {
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		ui->SessionStart();
		
		delete ui;
	}

	delete visManager;
	delete runManager;
}


