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
// $Id: MedLinac.cc,v 1.12 2007/07/01 15:19:18 mpiergen Exp $
//
// --------------------------------------------------------------
//      GEANT 4 -  medical_linac
//
// Code developed by: M. Piergentili

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4UImessenger.hh"
#include "Randomize.hh" 


#include "MedLinacDetectorMessenger.hh"
#include "MedLinacDetectorConstruction.hh"
#include "MedLinacPhysicsList.hh"
#include "MedLinacPhantomSD.hh"
#include "MedLinacPhantomHit.hh"
#include "MedLinacPrimaryGeneratorAction.hh"
#include "MedLinacEventAction.hh"
#include "MedLinacRunAction.hh"
#include "MedLinacTrackingAction.hh"
#include "G4SDManager.hh"


#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

int main(int argc ,char ** argv)
{

//choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  G4int seed = 0;
  CLHEP::HepRandom :: setTheSeed(seed);

 // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;
  G4String sensitiveDetectorName = "Phantom";
  MedLinacDetectorConstruction* pDetectorConstruction = 
MedLinacDetectorConstruction::GetInstance(sensitiveDetectorName);


  // set mandatory initialization classes
  runManager->SetUserInitialization(pDetectorConstruction);
  runManager->SetUserInitialization(new MedLinacPhysicsList);


#ifdef G4VIS_USE
  // Visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif


  // output environment variables:
#ifdef G4ANALYSIS_USE
  G4cout << G4endl << G4endl << G4endl 
	 << " User Environment " << G4endl
	 << " Using AIDA 3.0 analysis " << G4endl;
# else
  G4cout << G4endl << G4endl << G4endl 
	 << " User Environment " << G4endl
	 << " G4ANALYSIS_USE environment variable not set, NO ANALYSIS " 
	 << G4endl;
#endif
  



 G4UIsession* session = 0;
  if (argc == 1)   // Define UI session for interactive mode.
    {
      session = new G4UIterminal();
    }

  // set mandatory user action class
  
  runManager->SetUserAction(new MedLinacPrimaryGeneratorAction);
  MedLinacEventAction *pMedLinacEventAction = new MedLinacEventAction();
  runManager->SetUserAction(pMedLinacEventAction);

  MedLinacRunAction *pRunAction = new MedLinacRunAction();

  runManager->SetUserAction(pRunAction);  
  runManager->SetUserAction(new MedLinacTrackingAction);


  // Initialize  G4 kernel
  //runManager->Initialize();

  // get the pointer to the UI manager
  G4UImanager* UI = G4UImanager::GetUIpointer();
  
 if (session)   // Define UI session for interactive mode.
    { 
      G4cout<<" UI session starts ..."<< G4endl; 
       UI -> ApplyCommand("/control/execute macro.mac");    
       session->SessionStart();
       delete session;
    }

  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }  


  //job termination

#ifdef G4VIS_USE
 delete visManager;
#endif

 delete runManager;

 return 0;
}



