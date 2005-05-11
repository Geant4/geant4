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
// $Id: PhotIn.cc,v 1.1 2005-05-11 10:37:19 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - example/advanced/Photon_inefficiency 
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
#include "PhotInVisManager.hh"
#endif
#include "PhotInDetectorConstruction.hh"
#include "PhotInPhysicsList.hh"
#include "PhotInPrimaryGeneratorAction.hh"
#include "PhotInRunAction.hh"
#include "PhotInEventAction.hh"
#include "PhotInStackingAction.hh"

int main(int argc,char** argv) {

 // Construct the default run manager
 G4RunManager * runManager = new G4RunManager;

 // set mandatory initialization classes
 runManager->SetUserInitialization(new PhotInDetectorConstruction);
 runManager->SetUserInitialization(new PhotInPhysicsList);
  
#ifdef G4VIS_USE
 // visualization manager
 G4VisManager* visManager = new PhotInVisManager;
 visManager->Initialize();
#endif
    
 // set user action classes
 runManager->SetUserAction(new PhotInPrimaryGeneratorAction);
 runManager->SetUserAction(new PhotInRunAction);
 runManager->SetUserAction(new PhotInEventAction);
 runManager->SetUserAction(new PhotInStackingAction);
  
 // Initialize G4 kernel
 // runManager->Initialize();
    
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

