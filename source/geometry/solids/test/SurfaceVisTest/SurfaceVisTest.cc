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
//geotest.cc
//Geant4 includes
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4Circle.hh"
#include "G4VisExecutive.hh"
#include "G4VSolid.hh"
#include "globals.hh"
#include <vector>
#include "G4LogicalVolumeStore.hh"
// random number generator
#include "Randomize.hh" 

//a pre-built physics list
#include "QGSP_BERT.hh"

//My includes
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"


int main(int argc,char** argv)
{
        //Random number generator=always starts with the same seed
        CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); 
        CLHEP::HepRandom::setTheEngine( &defaultEngine ); 
        G4int seed = time( NULL ); 
        CLHEP::HepRandom::setTheSeed( seed ); 
        G4cout << G4endl 
         << " ===================================================== " << G4endl 
         << " Initial seed = " << seed << G4endl 
	 << " ===================================================== " << G4endl 
	 << G4endl; 

	//Construct the default run manager
	G4RunManager* runManager = new G4RunManager;
        G4VisManager* visManager = new G4VisExecutive;


  	//Set mandatory initialization and user action classes
 	runManager->SetUserInitialization(new QGSP_BERT);
       
  	runManager->SetUserInitialization(new DetectorConstruction);
  	runManager->SetUserAction(new PrimaryGeneratorAction);
        RunAction *runAction = new RunAction;
        runManager->SetUserAction(runAction);
        EventAction* eventAction = new EventAction(runAction);
        runManager->SetUserAction(eventAction);
  	//Initialize G4 kernel
  	runManager->Initialize();
        visManager->Initialize();

  	// run initialisation macro

      if ( argc==1 ) {   // Define UI session for interactive mode. 
	//Open a tcsh session: will stay there until the user types "exit"

  	G4UIsession* session = new G4UIterminal(new G4UItcsh);
         G4UImanager* UI = G4UImanager::GetUIpointer(); 
        UI->ApplyCommand("/control/execute vis.mac");
      
	
	session->SessionStart();
	delete session;
      }else{// Batch mode 
        G4String command = "/control/execute "; 
        G4String fileName = argv[1]; 
       	G4UIsession* session = new G4UIterminal(new G4UItcsh);
        G4UImanager* UI = G4UImanager::GetUIpointer(); 
        UI->ApplyCommand(command+fileName); 
        session->SessionStart();
	delete session;
      }
	
	delete visManager;
	
	// job termination
  	delete runManager;

  	return 0;

}// end of main()

//end of geotest.cc

