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
// $Id: ParN04.cc,v 1.5 2006/06/29 17:35:16 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - ParN04
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
#include "G4VisExecutive.hh"
#endif

#include "ParTopC.icc"

int main(int argc,char** argv)
{
  G4VSteppingVerbose::SetInstance(new ExN04SteppingVerbose);
  
  G4RunManager* runManager = new G4RunManager;

  runManager->SetUserInitialization(new ExN04DetectorConstruction);
  runManager->SetUserInitialization(new ExN04PhysicsList);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;
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

