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
// $Id: exampleRE02.cc,v 1.5 2010-11-08 18:48:54 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
//
#include "RE02DetectorConstruction.hh"
#include "RE02PhysicsList.hh"
#include "RE02PrimaryGeneratorAction.hh"
#include "RE02RunAction.hh"
#include "RE02EventAction.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//
int main(int argc,char** argv) {

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes (mandatory)
  //---
  //  Create Detector
  RE02DetectorConstruction* detector = new RE02DetectorConstruction;
  detector->SetPhantomSize(G4ThreeVector(200*mm,200*mm,400*mm)); //Default
  detector->SetNumberOfSegmentsInPhantom(100,100,200); //Default
  //  If your machine does not have enough memory,
  //  please try to reduce Number of Segements in phantom.
  //detector->SetNumberOfSegmentsInPhantom(1,1,100);  // For small memory size.
  //
  detector->SetLeadSegment(TRUE); // Default (Water and Lead)
  // Water and Lead segments are placed alternately, by defult.,
  //  If you want to simulation homogeneous water phantom, please set it FALSE.
  //detector->SetLeadSegment(FALSE); // Homogeneous water phantom
  //
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new RE02PhysicsList);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
   
  // UserAction classes
  runManager->SetUserAction(new RE02PrimaryGeneratorAction);
  runManager->SetUserAction(new RE02RunAction);  
  runManager->SetUserAction(new RE02EventAction);

  //Initialize G4 kernel
  runManager->Initialize();
      
  //get the pointer to the User Interface manager 
  G4UImanager * UI = G4UImanager::GetUIpointer();  

  if(argc==1)
  // Define (G)UI terminal for interactive mode  
  { 
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    UI->ApplyCommand("/control/execute vis.mac");    
    ui->SessionStart();
    delete ui;
  }
  else
  // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//

