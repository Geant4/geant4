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
// $Id: hepmcEx01.cc,v 1.2 2002-11-19 10:23:46 murakami Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN04
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------


#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "ExN04DetectorConstruction.hh"
#include "ExN04PhysicsList.hh"
#include "ExN04PrimaryGeneratorAction.hh"
#include "ExN04RunAction.hh"
#include "ExN04EventAction.hh"
#include "ExN04StackingAction.hh"
#include "ExN04TrackingAction.hh"
#include "ExN04SteppingAction.hh"
#include "ExN04SteppingVerbose.hh"

#ifdef G4VIS_USE
#include "ExN04VisManager.hh"
#endif

int main(int argc,char** argv)
{
  G4VSteppingVerbose::SetInstance(new ExN04SteppingVerbose);
  
  G4RunManager* runManager = new G4RunManager;

  runManager->SetUserInitialization(new ExN04DetectorConstruction);
  runManager->SetUserInitialization(new ExN04PhysicsList);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new ExN04VisManager;
  visManager->Initialize();
#endif

  runManager->Initialize();

  runManager->SetUserAction(new ExN04PrimaryGeneratorAction);
  runManager->SetUserAction(new ExN04RunAction);  
  runManager->SetUserAction(new ExN04EventAction);
  runManager->SetUserAction(new ExN04StackingAction);
  runManager->SetUserAction(new ExN04TrackingAction);
  runManager->SetUserAction(new ExN04SteppingAction);
  
  //get the pointer to the User Interface manager   
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if(argc==1)
  {
    // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_TCSH
    G4UIsession* session = new G4UIterminal(new G4UItcsh);      
#else
    G4UIsession* session = new G4UIterminal();
#endif    
    UImanager->ApplyCommand("/control/execute vis.mac");
    session->SessionStart();
    delete session;
  }
  else
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  delete runManager;
#ifdef G4VIS_USE
  delete visManager;
#endif

  return 0;
}

