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
// $Id: hepmcEx00.cc,v 1.1 2002-04-29 20:43:21 asaim Exp $
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
#include "ExN04EventAction.hh"
#include "ExN04StackingAction.hh"
#include "ExN04TrackingAction.hh"
#include "ExN04SteppingAction.hh"

#ifdef G4VIS_USE
#include "XXVisManager.hh"
#endif

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

#ifdef G4VIS_USE
  // initialize visualization package
  G4VisManager* visManager= new XXVisManager;
  visManager-> Initialize();
  G4cout << G4endl;
#endif

  if(argc==1) {
    G4UItcsh* tcsh= new
      G4UItcsh("[40;01;36mhepmcEx00[40;33m(%s)[40;32m[%/][00;30m:");
    G4UIterminal* session= new G4UIterminal(tcsh);
    tcsh-> SetLsColor(RED, GREEN);

    session->SessionStart();
    delete session;
  } else {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;

  return 0;
}

