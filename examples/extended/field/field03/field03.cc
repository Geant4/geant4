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
// $Id: field03.cc,v 1.4 2001-11-07 16:36:33 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - Example F03
//
// --------------------------------------------------------------
// Comments
//     
//   
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "F03VisManager.hh"
#endif

#include "F03DetectorConstruction.hh"
#include "F03ElectroMagneticField.hh"
#include "F03PhysicsList.hh"
#include "F03PrimaryGeneratorAction.hh"
#include "F03RunAction.hh"
#include "F03EventAction.hh"
#include "F03SteppingAction.hh"
#include "F03SteppingVerbose.hh"

int main(int argc,char** argv) 
{

  //choose the Random engine

  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class

  G4VSteppingVerbose::SetInstance(new F03SteppingVerbose);
  
  // Construct the default run manager

  G4RunManager * runManager = new G4RunManager;

  // Construct the electromagnetic field

  F03ElectroMagneticField* field = new F03ElectroMagneticField() ;
    
  // Set mandatory initialization classes

  F03DetectorConstruction* detector;
  detector = new F03DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new F03PhysicsList(detector));
  
#ifdef G4VIS_USE

  // visualization manager

  G4VisManager* visManager = new F03VisManager;
  visManager->Initialize();

#endif 
 
  // Set user action classes

  runManager->SetUserAction(new F03PrimaryGeneratorAction(detector));

  F03RunAction* runAction = new F03RunAction;

  runManager->SetUserAction(runAction);

  F03EventAction* eventAction = new F03EventAction(runAction);

  runManager->SetUserAction(eventAction);

  F03SteppingAction* steppingAction = new F03SteppingAction();
  runManager->SetUserAction(steppingAction);
  
  // Initialize G4 kernel, physics tables ...

  runManager->Initialize();
    
  // Get the pointer to the User Interface manager 

  G4UImanager* UI = G4UImanager::GetUIpointer();  
 
  if (argc==1)   // Define UI terminal for interactive mode  
  { 
     G4UIsession * session = new G4UIterminal;
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
  delete field;
  delete runManager;

  return 0;
}

