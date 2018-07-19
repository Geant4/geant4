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
// Please cite the following paper if you use this software
// Nucl.Instrum.Meth.B260:20-27, 2007

#ifdef G4MULTITHREADED
  #include "G4MTRunManager.hh"
#else
  #include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

int main(int argc,char** argv) {

  // Choose the Random engine
  
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager

#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;

  // Use only one thread for aberration coefficient calculation ("coef*" macros)
  //
  // For high statistics (no aberration coefficient calculation, "image*" & "grid*" macros), 
  // switch to more threads 
    
  runManager->SetNumberOfThreads(1);

#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  
  DetectorConstruction* detector = new DetectorConstruction;
    
  runManager->SetUserInitialization(detector);
  
  runManager->SetUserInitialization(new PhysicsList);
  
  // User action initialization
  
  runManager->SetUserInitialization(new ActionInitialization(detector));
    
  // Initialize G4 kernel
  
  runManager->Initialize();
    
  // Get the pointer to the User Interface manager 
  
  G4UImanager* UImanager = G4UImanager::GetUIpointer(); 
  
  if (argc==1)   // Define UI session for interactive mode.
  { 
#ifdef _WIN32
    G4UIsession * session = new G4UIterminal();
#else
    G4UIsession * session = new G4UIterminal(new G4UItcsh);
#endif
    UImanager->ApplyCommand("/control/execute default.mac");    
    session->SessionStart();
    delete session;
  }
  else           // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  //
  
  delete runManager;

  return 0;
}
