//
// fred: GEANT4 test program
//

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "FredDetectorConstruction.hh"
#include "FredPhysicsList.hh"
#include "FredPrimaryGeneratorAction.hh"
#include "FredVisManager.hh"
#include "FredEventAction.hh"
#include "FredMessenger.hh"

int main()
{
	// Construct run manager
	G4RunManager *runManager = new G4RunManager;
	
	// Build our master control messenger
	FredMessenger *messenger = new FredMessenger;
	
	// Build our detector
	runManager->SetUserInitialization( new FredDetectorConstruction(messenger) );
	runManager->SetUserInitialization( new FredPhysicsList );
	

	// Initialize visualization manager
	FredVisManager *visManager = new FredVisManager;
	visManager->Initialize();

	// Define our generator
	runManager->SetUserAction( new FredPrimaryGeneratorAction(messenger) );
	
	// Do something interesting at the end of each event
	runManager->SetUserAction( new FredEventAction(messenger) );
	
	// Give control to interactive terminal
	G4UIsession *session = new G4UIterminal;
	session->SessionStart();
	
	// All finished...
	delete session;
	delete runManager;
	delete visManager;
	return 0;
}
