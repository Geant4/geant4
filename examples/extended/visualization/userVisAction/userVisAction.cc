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
// $Id: userVisAction.cc,v 1.3 2009-12-22 15:07:19 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "UVA_DetectorConstruction.hh"
#include "UVA_PhysicsList.hh"
#include "UVA_PrimaryGeneratorAction.hh"
#include "UVA_RunAction.hh"
#include "UVA_EventAction.hh"
#include "UVA_SteppingAction.hh"
#include "UVA_SteppingVerbose.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#include "UVA_VisAction.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new UVA_SteppingVerbose);
  
  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes (mandatory)
  UVA_DetectorConstruction* UVA_detector = new UVA_DetectorConstruction;
  runManager->SetUserInitialization(UVA_detector);
  runManager->SetUserInitialization(new UVA_PhysicsList);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  // Set User Vis Action...
  visManager->SetUserAction
    (new UVA_VisAction, G4VisExtent(-1*m,1*m,-2*m,0,3.5*m,5.5*m));
  // 2nd argument optional - overridden by /vis/scene/add/userAction
  // arguments, if any.
#endif
   
  // UserAction classes
  runManager->SetUserAction(new UVA_PrimaryGeneratorAction(UVA_detector));
  runManager->SetUserAction(new UVA_RunAction);  
  runManager->SetUserAction(new UVA_EventAction);
  runManager->SetUserAction(new UVA_SteppingAction);

  //Initialize G4 kernel
  runManager->Initialize();
      
  //get the pointer to the User Interface manager 
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  

  if(argc==1)
#ifdef G4UI_USE
  // Define (G)UI terminal for interactive mode  
  { 
    G4UIExecutive * ui = new G4UIExecutive(argc, argv);
    UImanager->ApplyCommand("/control/execute vis.mac");    
    ui->SessionStart();
    delete ui;
#endif
  }
  else
  // Batch mode
  { 
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

