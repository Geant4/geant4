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
// $Id: SurfaceChecker.cc,v 1.1 2005-05-19 13:07:29 link Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "Randomize.hh"

#include "SCDetectorConstruction.hh"
#include "SCPhysicsList.hh"
#include "SCPrimaryGeneratorAction.hh"
#include "SCRunAction.hh"
#include "SCEventAction.hh"
#include "SCSteppingAction.hh"
#include "SCSteppingVerbose.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "SCVisManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  // choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);

  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SCSteppingVerbose);
  
  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes (mandatory)
  SCDetectorConstruction* SCdetector = new SCDetectorConstruction;
  runManager->SetUserInitialization(SCdetector);
  runManager->SetUserInitialization(new SCPhysicsList);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new SCVisManager;
  visManager->Initialize();
#endif
   
  // UserAction classes
  runManager->SetUserAction(new SCPrimaryGeneratorAction(SCdetector));
  runManager->SetUserAction(new SCRunAction);  
  runManager->SetUserAction(new SCEventAction);
  runManager->SetUserAction(new SCSteppingAction);

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

