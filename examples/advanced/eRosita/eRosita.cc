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
// $Id: eRosita.cc,v 1.2 2010-11-23 20:09:32 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "eRositaDetectorConstruction.hh"
#include "eRositaPhysicsList.hh"
#include "eRositaPrimaryGeneratorAction.hh"
#include "eRositaRunAction.hh"
#include "eRositaEventAction.hh"
#include "eRositaSteppingAction.hh"
#include "eRositaSteppingVerbose.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif


int main(int argc,char** argv)
{
  // User Verbose output class
  //
  G4VSteppingVerbose* verbosity = new eRositaSteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
  
  // Run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // User Initialization classes (mandatory)
  //
  eRositaDetectorConstruction* detector = new eRositaDetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new eRositaPhysicsList;
  runManager->SetUserInitialization(physics);
   
  // User Action classes
  //
  G4VUserPrimaryGeneratorAction* genAction = new eRositaPrimaryGeneratorAction(detector);
  runManager->SetUserAction(genAction);
  //
  G4UserRunAction* runAction = new eRositaRunAction;
  runManager->SetUserAction(runAction);
  //
  G4UserEventAction* eventAction = new eRositaEventAction;
  runManager->SetUserAction(eventAction);
  //
  G4UserSteppingAction* steppingAction = new eRositaSteppingAction;
  runManager->SetUserAction(steppingAction);

  // Initialize G4 kernel
  //
  runManager->Initialize();
      
  // Get the pointer to the User Interface manager
  //
  G4UImanager * UI = G4UImanager::GetUIpointer();  

  if (argc!=1)   // batch mode  
    {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
    }
    
  else           // interactive mode : define visualization and UI terminal
    { 
#ifdef G4VIS_USE
      G4VisManager* visManager = new G4VisExecutive;
      visManager->Initialize();
#endif    
     
      G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif
      UI->ApplyCommand("/control/execute vis.mac");     
      session->SessionStart();
      delete session;
     
#ifdef G4VIS_USE
      delete visManager;
#endif     
    }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  delete runManager;
  delete verbosity;

  return 0;
}


