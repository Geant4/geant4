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
// $Id: bench02.cc,v 1.2 2010-11-30 11:50:09 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "Em2DetectorConstruction.hh"
#include "Em2PhysicsList.hh"
#include "Em2PrimaryGeneratorAction.hh"
#include "Em2RunAction.hh"
#include "Em2EventAction.hh"
#include "Em2TrackingAction.hh"
#include "Em2SteppingAction.hh"
#include "Em2SteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {
 
  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new Em2SteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Em2DetectorConstruction* detector = new Em2DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Em2PhysicsList);
  
  Em2PrimaryGeneratorAction* primary = new Em2PrimaryGeneratorAction(detector);
  runManager->SetUserAction(primary);
    
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new Em2VisManager;
  visManager->Initialize();
#endif
    
  // set user action classes
  Em2RunAction* RunAct = new Em2RunAction(detector,primary);
  runManager->SetUserAction(RunAct);
  runManager->SetUserAction(new Em2EventAction   (RunAct));
  runManager->SetUserAction(new Em2TrackingAction(RunAct));
  runManager->SetUserAction(new Em2SteppingAction(detector,RunAct)); 
  
  // get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if (argc==1)   // Define UI terminal for interactive mode.
    {
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      ui->SessionStart();
      delete ui;
#endif
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

