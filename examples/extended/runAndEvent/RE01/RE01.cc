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
// $Id: RE01.cc,v 1.2 2005/06/01 18:27:29 perl Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// --------------------------------------------------------------
//      GEANT4 - RE01 exsample code
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

#include "RE01DetectorConstruction.hh"
#include "RE01PhysicsList.hh"
#include "RE01PrimaryGeneratorAction.hh"
#include "RE01RunAction.hh"
#include "RE01EventAction.hh"
#include "RE01StackingAction.hh"
#include "RE01TrackingAction.hh"
#include "RE01SteppingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

int main(int argc,char** argv)
{
  G4RunManager* runManager = new G4RunManager;

#ifdef G4VIS_USE
  // Visualization manager construction
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  runManager->SetUserInitialization(new RE01DetectorConstruction);
  runManager->SetUserInitialization(new RE01PhysicsList);
  
  runManager->Initialize();

  runManager->SetUserAction(new RE01PrimaryGeneratorAction);
  runManager->SetUserAction(new RE01RunAction);  
  runManager->SetUserAction(new RE01EventAction);
  runManager->SetUserAction(new RE01StackingAction);
  runManager->SetUserAction(new RE01TrackingAction);
  runManager->SetUserAction(new RE01SteppingAction);
  
  if(argc==1)
  {
#ifdef G4UI_USE_TCSH
    G4UIsession* session = new G4UIterminal(new G4UItcsh);      
#else
    G4UIsession* session = new G4UIterminal();
#endif    
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

#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;

  return 0;
}

