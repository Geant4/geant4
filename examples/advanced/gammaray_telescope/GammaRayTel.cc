// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTel.cc,v 1.3 2000-12-06 16:53:12 flongo Exp $
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
//      ------------ GammaRayTel example main program ------
//           by F.Longo, R.Giannitrapani & G.Santin (29 nov 2000)
//           See README file for details on this example            
// ************************************************************

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4VIS_USE
#include "GammaRayTelVisManager.hh"
#endif

#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelPhysicsList.hh"
#include "GammaRayTelPrimaryGeneratorAction.hh"
#include "GammaRayTelRunAction.hh"
#include "GammaRayTelEventAction.hh"

#ifdef G4ANALYSIS_USE
#include "GammaRayTelAnalysisManager.hh"
#endif


/* This global file is used to store relevant data for
   analysis with external tools */
G4std::ofstream outFile;

// This is the main function 
int main(int argc, char** argv)
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;
  
  // Set mandatory user initialization classes
  GammaRayTelDetectorConstruction* detector = 
    new GammaRayTelDetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new GammaRayTelPhysicsList);

  // Set mandatory user action classes
  runManager->SetUserAction(new GammaRayTelPrimaryGeneratorAction(detector));


#ifdef G4ANALYSIS_USE
  // Creation of the analysis manager
  GammaRayTelAnalysisManager* analysisMgr = new GammaRayTelAnalysisManager(detector);
#endif

  // Set optional user action classes
#ifdef G4ANALYSIS_USE
  GammaRayTelEventAction* eventAction = 
    new GammaRayTelEventAction(analysisMgr);
  GammaRayTelRunAction* runAction =
    new GammaRayTelRunAction(analysisMgr);
#else 
  GammaRayTelEventAction* eventAction = new GammaRayTelEventAction();
  GammaRayTelRunAction* runAction = new GammaRayTelRunAction();
#endif
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
#ifdef G4VIS_USE
  // Visualization manager
  G4VisManager* visManager = new GammaRayTelVisManager;
  visManager->Initialize();
#endif
  
  // Initialize G4 kernel
  runManager->Initialize();
  
  // Get the pointer to the UI manager
  G4UImanager* UI = G4UImanager::GetUIpointer();
  
  if (session) 
    {
      /* prerunGammaRayTel.mac is loaded by default
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
      else  UI->ApplyCommand("/control/execute prerunGammaRayTel.mac");
      session->SessionStart();
      delete session;
    }

  // Job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
#ifdef G4ANALYSIS_USE
  delete analysisMgr;
#endif
  delete runManager;
  return 0;
}








