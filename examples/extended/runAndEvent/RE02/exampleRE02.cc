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
// $Id: exampleRE02.cc,v 1.1 2005/11/24 01:44:18 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

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
    //  If your machine does not have enough memory,
    //  please try to reduce Number of Segements in phantom.
  RE02DetectorConstruction* detector = new RE02DetectorConstruction;
  detector->SetPhantomSize(G4ThreeVector(200*mm,200*mm,400*mm)); //Default
  detector->SetNumberOfSegmentsInPhantom(100,100,200); //Default
  //detector->SetNumberOfSegmentsInPhantom(1,1,100);  // For small memory size.
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
    // G4UIterminal is a (dumb) terminal.
    G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif    

    UI->ApplyCommand("/control/execute vis.mac");    
    session->SessionStart();
    delete session;
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

