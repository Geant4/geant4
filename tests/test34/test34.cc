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

#include <algorithm>
#include <iostream>
#include "G4Timer.hh"

// G4 includes 
#include "G4ios.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "G4Timer.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

// my project 
#include "Tst34DetectorConstruction.hh"
#include "Tst34PhysicsList.hh"
#include "Tst34PrimaryGeneratorAction.hh"
#include "Tst34EventAction.hh"
#include "Tst34RunAction.hh"

G4Timer Timer;
G4Timer Timerintern;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Timer to see GFlash performance
  Timer.Start();

  G4cout<<"+-------------------------------------------------------+"<<G4endl;
  G4cout<<"|                                                       |"<<G4endl;
  G4cout<<"|          This is an example of Shower                 |"<<G4endl;
  G4cout<<"|          Parameterization with GFLASH                 |"<<G4endl;
  G4cout<<"+-------------------------------------------------------+"<<G4endl;

  G4RunManager* runManager = new G4RunManager;

  // UserInitialization classes (mandatory)
  G4cout<<"# GFlash Example: Detector Construction"<<G4endl;    
  runManager->SetUserInitialization(new Tst34DetectorConstruction);
  G4cout<<"# GFlash Example: Physics list"<<G4endl;
  runManager->SetUserInitialization(new Tst34PhysicsList);
  G4cout<<"# GFlash Example: Primary Generator"<<G4endl;
  runManager->SetUserAction(new Tst34PrimaryGeneratorAction);
  G4cout<<"# GFlash Example: User Action Classes"<<G4endl;
  runManager->SetUserAction(new Tst34EventAction);
  runManager->SetUserAction(new Tst34RunAction);

  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/run/verbose 0");
  runManager->Initialize();
  UI->ApplyCommand("/Step/Verbose 0");

  if (argc==1)   // Define UI terminal for interactive mode  
  { 
    // G4UIterminal is a (dumb) terminal.
    G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif    
    session->SessionStart();
    delete session;
  }
  else           // Batch mode
  { 
    G4String s=*(argv+1);
    UI->ApplyCommand("/control/execute "+s);
  }

  delete runManager;

  Timer.Stop();
  G4cout << G4endl;
  G4cout << "******************************************";
  G4cout << G4endl;
  G4cout << "Total Real Elapsed Time is: "<< Timer.GetRealElapsed();
  G4cout << G4endl;
  G4cout << "Total System Elapsed Time: " << Timer.GetSystemElapsed();
  G4cout << G4endl;
  G4cout << "Total GetUserElapsed Time: " << Timer.GetUserElapsed();
  G4cout << G4endl;
  G4cout << "******************************************";
  G4cout << G4endl;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
