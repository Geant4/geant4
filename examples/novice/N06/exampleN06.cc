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
// $Id: exampleN06.cc,v 1.7 2003-01-23 15:31:39 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Description: Test of Continuous Process G4Cerenkov
//              and RestDiscrete Process G4Scintillation
//              -- Generation Cerenkov Photons --
//              -- Generation Scintillation Photons --
//              -- Transport of optical Photons --
// Version:     5.0
// Created:     1996-04-30
// Author:      Juliet Armstrong
// mail:        gum@triumf.ca

// Updated:
//     
// 2002-11-12 by Peter Gumplinger
//   - Change user interface to G4Scintillation
//   - Add UserStackingAction
//
// 1998-06-03 by Peter Gumplinger for example ExN06
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "G4ios.hh"
#include <stdlib.h>

#include "ExN06DetectorConstruction.hh"
#include "ExN06PhysicsList.hh"
#include "ExN06PrimaryGeneratorAction.hh"
#include "ExN06RunAction.hh"
#include "ExN06EventAction.hh"
#include "ExN06StackingAction.hh"
#include "ExN06SteppingVerbose.hh"

#ifdef G4VIS_USE
#include "ExN06VisManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  // Seed the random number generator manually
  G4long myseed = 345354;
  HepRandom::setTheSeed(myseed);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new ExN06SteppingVerbose);
  
  // Run manager
  G4RunManager* runManager = new G4RunManager;

  // UserInitialization classes - mandatory
  runManager-> SetUserInitialization(new ExN06DetectorConstruction);
  runManager-> SetUserInitialization(new ExN06PhysicsList);
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new ExN06VisManager;
  visManager->Initialize();
#endif

  // UserAction classes 
  runManager->SetUserAction(new ExN06RunAction);
  runManager->SetUserAction(new ExN06PrimaryGeneratorAction);
  runManager->SetUserAction(new ExN06EventAction);  
  runManager->SetUserAction(new ExN06StackingAction);
  
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer(); 
   
  if (argc==1)   //define UI session for interactive mode
    {
      G4UIsession* session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif    
      UI->ApplyCommand("/control/execute vis.mac"); 
      session->SessionStart();
      delete session;
   }
   
  else         // Batch mode
   {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
   }
   
  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
