// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test20.cc,v 1.1 2001-05-24 19:48:58 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 main program
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Test20 main program ------
//           by G.Depaola & F.Longo (may 2001) 
// ************************************************************

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#include "Tst20DetectorConstruction.hh"
#include "Tst20PhysicsList.hh"
#include "Tst20PrimaryGeneratorAction.hh"
#include "Tst20RunAction.hh"
#include "Tst20EventAction.hh"

// This is the main function 

int main(int argc, char** argv)
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;
  
  // Set mandatory user initialization classes
  Tst20DetectorConstruction* detector = 
    new Tst20DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Tst20PhysicsList);

  // Set mandatory user action classes
  runManager->SetUserAction(new Tst20PrimaryGeneratorAction(detector));

  Tst20EventAction* eventAction = new Tst20EventAction();
  Tst20RunAction* runAction = new Tst20RunAction();
  
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(runAction);
  
  // Set visualization and user interface
  // Initialization of the User Interface Session

  G4UIsession* session=0;

#ifdef G4UI_USE_XM
  // Create a XMotif user interface
  session = new G4UIXm(argc,argv);
#else
  // Create the standard user interface
  session = new G4UIterminal;
#endif
  
  // Initialize G4 kernel
  runManager->Initialize();
  
  // Get the pointer to the UI manager
  G4UImanager* UI = G4UImanager::GetUIpointer();
  
  if (session) 
    {
      /* prerunTst20.mac is loaded by default
	 unless a macro file is passed as the argument
	 of the executable */

      if(argc>1)
	{
	  G4String command = "/control/execute ";
	  for (int i=2; i<=argc; i++) 
	    {
	      G4String macroFileName = argv[i-1];
	      UI->ApplyCommand(command+macroFileName);
	    }
	}
      else  UI->ApplyCommand("/control/execute prerunTst20.mac");
      session->SessionStart();
      delete session;
    }

  delete runManager;
  return 0;
}








