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
// $Id: GammaRayTel.cc,v 1.7 2001-11-29 11:19:16 griccard Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 main program
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTel example main program ------
//           by F.Longo, R.Giannitrapani & G.Santin (29 nov 2000)
//           See README file for details on this example            
//  20.11.01 G.Santin: new analysis management, and some modification in the 
//                     construction of some Action's
// ************************************************************

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

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
#include "GammaRayTelAnalysis.hh"
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
  runManager->SetUserAction(new GammaRayTelPrimaryGeneratorAction);


#ifdef G4ANALYSIS_USE
  // Creation of the analysis manager
  GammaRayTelAnalysis* analysis = GammaRayTelAnalysis::getInstance();
#endif

  // Set optional user action classes
  GammaRayTelEventAction* eventAction = new GammaRayTelEventAction();
  GammaRayTelRunAction* runAction = new GammaRayTelRunAction();
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(runAction);

  // Set visualization and user interface
  // Initialization of the User Interface Session
  G4UIsession* session=0;
#ifdef G4UI_USE_XM
  // Create a XMotif user interface
  session = new G4UIXm(argc,argv);
#else
#ifdef G4UI_USE_TCSH
  session = new G4UIterminal(new G4UItcsh);      
#else
  // Create the standard user interface
  session = new G4UIterminal();
#endif
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
  delete analysis;
#endif
  delete runManager;
  return 0;
}








