//geotest.cc
//Geant4 includes
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
#include "CLHEP/Random/RanluxEngine.h" 

//a pre-built physics list
#include "QGSP_EMV.hh"

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
 	runManager->SetUserInitialization(new QGSP_EMV);
       
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

