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
// $Id: A01app.cc,v 1.3 2002-12-13 11:34:27 gunter Exp $
// --------------------------------------------------------------
//
// --------------------------------------------------------------
//      GEANT 4 - A01app
// --------------------------------------------------------------
// Comments
//   Tutorial for Geant4 lectures
// --------------------------------------------------------------


#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "A01DetectorConstruction.hh"
#include "A01PhysicsList.hh"
#include "A01PrimaryGeneratorAction.hh"
#include "A01EventAction.hh"
#include "A01TrackingAction.hh"

#ifdef G4VIS_USE
#include "A01VisManager.hh"
#endif

int main(int argc,char** argv)
{
  // RunManager construction
  G4RunManager* runManager = new G4RunManager;

#ifdef G4VIS_USE
  // Visualization manager construction
  G4VisManager* visManager = new A01VisManager();
  visManager->Initialize();
#endif

  // mandatory user initialization classes
  runManager->SetUserInitialization(new A01DetectorConstruction);
  runManager->SetUserInitialization(new A01PhysicsList);

  // initialize Geant4 kernel
  runManager->Initialize();

  // mandatory user action class
  runManager->SetUserAction(new A01PrimaryGeneratorAction);

  // optional user action classes
  runManager->SetUserAction(new A01EventAction);
  runManager->SetUserAction(new A01TrackingAction);

  if(argc>1)
  // execute an argument macro file if exist
  {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else
  // start interactive session
  {
    G4UIsession* session = new G4UIterminal();
    session->SessionStart();
    delete session;
  }

#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;

  return 0;
}

