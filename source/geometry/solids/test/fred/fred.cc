//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// fred: GEANT4 test program
//

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "FredDetectorConstruction.hh"
#include "FredPhysicsList.hh"
#include "FredPrimaryGeneratorAction.hh"
#include "FredVisManager.hh"
#include "FredEventAction.hh"
#include "FredMessenger.hh"

int main(int argc,char *argv[])
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
  G4UIsession *session;
  
#ifdef G4UI_USE_TCSH
  session = new G4UIterminal(new G4UItcsh);      
#else
	session = new G4UIterminal;
#endif
  
  if (argc > 1)  {
    G4UImanager * UI = G4UImanager::GetUIpointer();

    for (int i=1;i<argc;i++)    {
      UI->ApplyCommand("/control/execute "+G4String(argv[i]));
    }
  }
  else  {
    session->SessionStart();
  }
	
	// All finished...
	delete session;
	delete runManager;
	delete visManager;
	return 0;
}
