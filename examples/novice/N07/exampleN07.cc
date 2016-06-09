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
// $Id: exampleN07.cc,v 1.3 2003/04/25 17:03:36 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN07 
//
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "ExN07VisManager.hh"
#endif
#include "ExN07DetectorConstruction.hh"
#include "ExN07PhysicsList.hh"
#include "ExN07PrimaryGeneratorAction.hh"
#include "ExN07RunAction.hh"
#include "ExN07EventAction.hh"
#include "ExN07StackingAction.hh"

int main(int argc,char** argv) {

 // Construct the default run manager
 G4RunManager * runManager = new G4RunManager;

 // set mandatory initialization classes
 runManager->SetUserInitialization(new ExN07DetectorConstruction);
 runManager->SetUserInitialization(new ExN07PhysicsList);
  
#ifdef G4VIS_USE
 // visualization manager
 G4VisManager* visManager = new ExN07VisManager;
 visManager->Initialize();
#endif
    
 // set user action classes
 runManager->SetUserAction(new ExN07PrimaryGeneratorAction);
 runManager->SetUserAction(new ExN07RunAction);
 runManager->SetUserAction(new ExN07EventAction);
 runManager->SetUserAction(new ExN07StackingAction);
  
 //Initialize G4 kernel
 runManager->Initialize();
    
 // get the pointer to the User Interface manager 
 G4UImanager* UI = G4UImanager::GetUIpointer();  

 if (argc==1)   // Define UI session for interactive mode.
 {
   G4UIsession* session=0;
  
   // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_TCSH
   session = new G4UIterminal(new G4UItcsh);      
#else
   session = new G4UIterminal();
#endif    
      
   UI->ApplyCommand("/control/execute vis.mac");    
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
 delete runManager;

 return 0;
}

