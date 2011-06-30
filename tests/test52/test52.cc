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

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGenerator.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "globals.hh"
// #ifdef G4VIS_USE
// #include "G4VisExecutive.hh"
// #include "G4TrajectoryGenericDrawer.hh"
// #endif


int main(int argc,char** argv) {
 	
  G4RunManager * runManager = new G4RunManager;

  DetectorConstruction* detector = new DetectorConstruction();
  runManager -> SetUserInitialization(detector);  

  PhysicsList* physics = new PhysicsList();
  runManager -> SetUserInitialization(physics);

  PrimaryGenerator* source = new PrimaryGenerator(); 
  runManager -> SetUserAction(source);

  RunAction* runAction = new RunAction(); 
  runManager -> SetUserAction(runAction);  

  EventAction* eventAction = new EventAction();
  runManager -> SetUserAction(eventAction);
     
  SteppingAction* steppingAction = new SteppingAction();
  runManager -> SetUserAction(steppingAction);

  TrackingAction* trackingAction = new TrackingAction();
  runManager -> SetUserAction(trackingAction);

// #ifdef G4VIS_USE 
//   std::cout << "INFORMATION: Visualisation used." << std::endl; 

//   G4VisManager* visManager = new G4VisExecutive;
//   visManager -> Initialize();

//   G4TrajectoryGenericDrawer* genDrawer = new G4TrajectoryGenericDrawer;
//   visManager -> RegisterModel(genDrawer);
// #endif

  G4UImanager * UI = G4UImanager::GetUIpointer();  
   
  if(argc != 2) {
     std::cerr << "Error. Wrong number of arguments." << std::endl;
     std::cout << "INFORMATION: Usage " << argv[0] 
               << " macrofile" << std::endl; 
  }
  else if(argc == 2) { 
     G4String fileName = argv[1];
     std::cout << "INFORMATION:  Commands in file " << fileName 
               << " used to control simulation." << std::endl; 
     UI -> ApplyCommand("/control/execute "+fileName);
  }

// #ifdef G4VIS_USE
//   delete visManager;
// #endif

  delete runManager;
  return 0;
}
