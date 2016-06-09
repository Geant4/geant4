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
// $Id: lArCal.cc,v 1.13 2010-11-16 17:54:39 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN03 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
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

#include "FCALTestbeamSetup.hh"
#include "FCALSteppingVerbose.hh"
#include "FCALPrimaryGeneratorAction.hh"
#include "LHEP.hh"
#include "QGSP.hh"
#include "QGSP_BIC_EMY.hh"


#ifdef G4ANALYSIS_USE

#include "FCALRunAction.hh"
#include "FCALTBEventAction.hh"
#include "FCALSteppingAction.hh"

#endif


int main(int argc,char** argv) {

  // choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new FCALSteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  FCALTestbeamSetup* detector = new FCALTestbeamSetup;
  runManager->SetUserInitialization(detector);

  //***LOOKHERE*** CHOOSE THE PHYSICS LIST.
  // runManager->SetUserInitialization(new LHEP);     // LHEP     
  // runManager->SetUserInitialization(new QGSP);     // QGSP   
  // runManager->SetUserInitialization(new QGSC);     // QGSC
  runManager->SetUserInitialization(new QGSP_BIC_EMY);     // QGSP_BIC_EMY    
  //***endLOOKHERE***
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
    
  // set user action classes
  runManager->SetUserAction(new FCALPrimaryGeneratorAction());

#ifdef G4ANALYSIS_USE

  FCALRunAction* RunAction = new FCALRunAction;
  runManager ->SetUserAction(RunAction);

  FCALSteppingAction* StepAction = new FCALSteppingAction;
  runManager->SetUserAction(StepAction);

  //  runManager->SetUserAction(new FCALRunAction);
  
  runManager->SetUserAction(new FCALTBEventAction(StepAction));
  
  
#endif

  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if (argc==1)   // Define UI session for interactive mode.
     {
#ifdef G4UI_USE
       G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
       UImanager->ApplyCommand("/control/execute prerunlArcal.mac");     
#endif
       if (ui->IsGUI())
	 UImanager->ApplyCommand("/control/execute gui.mac");     
       ui->SessionStart();
       delete ui;
#endif
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

