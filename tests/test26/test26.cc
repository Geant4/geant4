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
// $Id: test26.cc,v 1.2 2003-02-01 18:14:57 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"


#include "Tst26DetectorConstruction.hh"
#include "Tst26PhysicsList.hh"
#include "Tst26PrimaryGeneratorAction.hh"
#include "Tst26RunAction.hh"
#include "Tst26EventAction.hh"
#include "Tst26TrackingAction.hh"
#include "Tst26SteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {
 
  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
       
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Tst26DetectorConstruction* detector = new Tst26DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Tst26PhysicsList);
  
  Tst26PrimaryGeneratorAction* primary = new Tst26PrimaryGeneratorAction(detector);
  runManager->SetUserAction(primary);
        
  // set user action classes
  Tst26RunAction* RunAct = new Tst26RunAction(detector,primary);
  runManager->SetUserAction(RunAct);
  runManager->SetUserAction(new Tst26EventAction   (RunAct));
  runManager->SetUserAction(new Tst26TrackingAction(RunAct));
  runManager->SetUserAction(new Tst26SteppingAction(detector,RunAct)); 
  
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc==1)   // Define UI terminal for interactive mode.
    {
      G4UIsession * session = 0;
      session = new G4UIterminal();
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

