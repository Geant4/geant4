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
// SBT (Solids Batch Test) main program
//

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4InteractiveSolid.hh"
#include "SBTMessenger.hh"
#include "SBTvoxelMessenger.hh"
#include "SBTVisManager.hh"

/*
MEDERNACH Emmanuel
Aug 2000

You could now run SBT with an argument script
and exit SBT.
*/

int main(int argc,char *argv[])
{
        // Initialize visualization manager
        SBTVisManager *visManager = new SBTVisManager;
        visManager->Initialize();

	// Build our "interactive" volume,
	// a volume that is also a messenger
	G4InteractiveSolid *interactiveSolid = new G4InteractiveSolid( "/solid/" );
	
	// Build our batch test messenger.
	// A new batch test will be created by it.
	// The test messenger gets the target solid
	// from a solid query class G4QuerySolid, as specified in the
	// second argument
	SBTMessenger *testMessenger = new SBTMessenger( "/test/", (G4SolidQuery *)interactiveSolid, visManager );

	//
	// Build our voxel test messenger
	SBTvoxelMessenger *voxelMessenger = new SBTvoxelMessenger( "/voxel/", (G4SolidQuery *)interactiveSolid, visManager );
		
	// Give control to interactive terminal
	G4UIsession *session = new G4UIterminal;
	   
	if (argc > 1)  // when run with an argument, run each scripts and exit
	  {
	    G4UImanager * UI = G4UImanager::GetUIpointer();

	    for (int i=1;i<argc;i++)
	      {
		UI->ApplyCommand("/control/execute "+G4String(argv[i]));
	      }
	  }
	else
	  {
	    session->SessionStart();
	  }
	
	// All finished...
	delete session;
	delete testMessenger;
	delete voxelMessenger;
	delete interactiveSolid;
	return 0;
}
