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
// $Id: exampleRMC01.cc,v 1.1 2009/11/19 22:41:18 ldesorgh Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleRMC01
//
// --------------------------------------------------------------
// Comments
//
// This example intends to show how to use the adjoint/reverse simulation mode.
// 
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4UI_USE_WIN32
#include "G4UIWin32.hh"
#endif

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif


#include "G4AdjointSimManager.hh"
#include "RMC01DetectorConstruction.hh"
#include "RMC01PrimaryGeneratorAction.hh"
#include "RMC01EventAction.hh"
#include "RMC01RunAction.hh"
#include "G4AdjointPhysicsList.hh"


#include "G4AdjointSimManager.hh"
#include "RMC01AdjointEventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {
   
   
 
  // Construct the default run manager
  G4RunManager * theRunManager = new G4RunManager;
   
  
  RMC01DetectorConstruction* detector = new RMC01DetectorConstruction();
  
  //Physics and geometry are declared as in a normal G4 application
  //--------------------------------------------
  theRunManager->SetUserInitialization(detector);
  theRunManager->SetUserInitialization(new G4AdjointPhysicsList); // in this physics list all the adjoint
  								  //  processes have to be declared
			
							  
  
  theRunManager->SetUserAction(new RMC01PrimaryGeneratorAction);
  theRunManager->SetUserAction(new RMC01EventAction);
  RMC01RunAction* theRunAction = new RMC01RunAction;
  theRunManager->SetUserAction(theRunAction);
 
  
  //The adjoint simulation manager will control the Reverse MC mode 
  //---------------------------------------------------------------
  
  G4AdjointSimManager* theAdjointSimManager = G4AdjointSimManager::GetInstance();
  
  //It is possible to define action that will be used duirng the adjoint tracking phase
  //
  
  theAdjointSimManager->SetAdjointRunAction(theRunAction);
  theAdjointSimManager->SetAdjointEventAction(new RMC01AdjointEventAction);
  
  G4UIsession* session=0;
  
  if (argc==1)   // Define UI session for interactive mode.
  {
     session = new G4UIterminal(new G4UItcsh);
	
  }
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  
#endif

  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (session)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.  
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete theRunManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
