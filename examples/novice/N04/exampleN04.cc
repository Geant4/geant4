// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: exampleN04.cc,v 1.1 1999-01-07 16:05:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN04
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------


#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "ExN04DetectorConstruction.hh"
#include "ExN04PhysicsList.hh"
#include "ExN04PrimaryGeneratorAction.hh"
#include "ExN04EventAction.hh"
#include "ExN04StackingAction.hh"
#include "ExN04TrackingAction.hh"
#include "ExN04SteppingAction.hh"


int main(int argc,char** argv)
{
  G4RunManager* runManager = new G4RunManager;

  runManager->SetUserInitialization(new ExN04DetectorConstruction);
  runManager->SetUserInitialization(new ExN04PhysicsList);

  runManager->Initialize();

  runManager->SetUserAction(new ExN04PrimaryGeneratorAction);
  runManager->SetUserAction(new ExN04EventAction);
  runManager->SetUserAction(new ExN04StackingAction);
  runManager->SetUserAction(new ExN04TrackingAction);
  runManager->SetUserAction(new ExN04SteppingAction);

  if(argc==1)
  {
    // G4UIterminal is a (dumb) terminal.
    G4UIsession* session = new G4UIterminal;
    session->SessionStart();
    delete session;
  }
  else
  {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  delete runManager;

  return 0;
}

