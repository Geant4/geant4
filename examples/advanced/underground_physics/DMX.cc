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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//            Underground Advanced example main program
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// main program
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "DMXAnalysisManager.hh"
#include "DMXDetectorConstruction.hh"
#include "DMXPhysicsList.hh"
#include "DMXPrimaryGeneratorAction.hh"
#include "DMXRunAction.hh"
#include "DMXEventAction.hh"
#include "DMXSteppingAction.hh"
#include "DMXStackingAction.hh"

#include <vector>

int main(int argc,char** argv) {

  // choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  DMXDetectorConstruction* detector = new DMXDetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new DMXPhysicsList);
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // output environment variables:
#ifdef G4ANALYSIS_USE
  G4cout << G4endl << G4endl << G4endl 
	 << " User Environment " << G4endl
	 << " Using AIDA 3.2.1 analysis " << G4endl;
#else
  G4cout << G4endl << G4endl << G4endl 
	 << " User Environment " << G4endl
	 << " G4ANALYSIS_USE environment variable not set, NO ANALYSIS " 
	 << G4endl;
#endif

#ifdef DMXENV_GPS_USE
  G4cout << " Using GPS and not DMX gun " << G4endl;
#else
  G4cout << " Using the DMX gun " << G4endl;
#endif
    
  // set user action classes
  DMXPrimaryGeneratorAction* DMXGenerator = new DMXPrimaryGeneratorAction;
  runManager->SetUserAction(DMXGenerator);
  //  runManager->SetUserAction(new DMXPrimaryGeneratorAction);
  // RunAction is inherited by EventAction for output filenames - will all
  // change when implement proper analysis manager?
  DMXRunAction* DMXRun = new DMXRunAction;
  runManager->SetUserAction(DMXRun);
  DMXEventAction* eventAction = new DMXEventAction(DMXRun,DMXGenerator);
  runManager->SetUserAction(eventAction);
  // eventAction is inherited by SteppingAction in order to switch colour
  // flag:
  runManager->SetUserAction(new DMXSteppingAction(eventAction));
  runManager->SetUserAction(new DMXStackingAction);

  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  // Define UI session for interactive mode.
  if(argc == 1)
    {
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute initInter.mac");     
#endif
      ui->SessionStart();
      delete ui;
#endif
    }
  // Batch mode
  else
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

  // job termination
#ifdef G4ANALYSIS_USE  
  DMXAnalysisManager::getInstance()->Finish(); 
  G4cout << "Analysis files closed" << G4endl;
#endif

#ifdef G4VIS_USE
  if(visManager) delete visManager;
#endif
  delete runManager;

  return 0;
}

