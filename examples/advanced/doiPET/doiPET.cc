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

//GEANT4 - Depth-of-Interaction enabled Positron emission tomography (PET) advanced example 

//Authors and contributors

// Author list to be updated, with names of co-authors and contributors from National Institute of Radiological Sciences (NIRS)

// Abdella M. Ahmed (1, 2), Andrew Chacon (1, 2), Harley Rutherford (1, 2),
// Hideaki Tashima (3), Go Akamatsu (3), Akram Mohammadi (3), Eiji Yoshida (3), Taiga Yamaya (3)
// Susanna Guatelli (2), and Mitra Safavi-Naeini (1, 2)

// (1) Australian Nuclear Science and Technology Organisation, Australia
// (2) University of Wollongong, Australia
// (3) National Institute of Radiological Sciences, Japan

#include "doiPETDetectorConstruction.hh"
#include "doiPETPhysicsList.hh"
#include "doiPETAnalysis.hh"
#include "doiPETActionInitialization.hh"
#include "Randomize.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManagerFactory.hh"
//
//////////////////////////////////////////////////////////////////////////////

int main(int argc,char** argv)
{
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  G4int nThreads = 4;
  runManager->SetNumberOfThreads(nThreads);
   
	G4UIExecutive* ui = 0;
	if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
	}
	
	//Initialize analysis
	doiPETAnalysis* ptrAnalysis = doiPETAnalysis::GetInstance();

	runManager->SetUserInitialization(new doiPETDetectorConstruction());

	runManager->SetUserInitialization(new doiPETPhysicsList());

	// Set user action initialization
	runManager->SetUserInitialization(new doiPETActionInitialization(ptrAnalysis));

	G4double act  = 1000000 * becquerel;//Activity is set via run.mac file
	ptrAnalysis->SetActivity(act);

	G4double halfLife = 109.771 * 60 * s; //Halflife of F-18 as a default
	ptrAnalysis -> SetIsotopeHalfLife(halfLife);

	//Blurring specification of the scanner. see inputParameter.txt
	ptrAnalysis -> BlurringParameters();

	//Open file to write the output of the simulation
	ptrAnalysis->Open("result"); //file extention is affixed based on the type of the output (.root for root or .data for ascii)


	//
	ptrAnalysis -> PMTPosition();
	//Read reflector pattern from the inputParameter.txt file
	ptrAnalysis->ReadReflectorPattern();

	// Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else { 
    // interactive mode
	G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
	delete visManager;
  }

	//close the file
	ptrAnalysis->Close();
	ptrAnalysis->Delete();
	
	delete runManager;
	return 0;

}

//////////////////////////////////////////////////////////////////////////////

