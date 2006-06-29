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
