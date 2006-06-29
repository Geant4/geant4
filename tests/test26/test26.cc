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
// $Id: test26.cc,v 1.7 2006-06-29 21:53:03 gunter Exp $
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
#include "G4UItcsh.hh"
#include "Randomize.hh"


#include "Tst26DetectorConstruction.hh"
#include "Tst26PhysicsList.hh"
#include "Tst26PrimaryGeneratorAction.hh"
#include "Tst26RunAction.hh"
#include "Tst26EventAction.hh"
#include "Tst26TrackingAction.hh"
#include "Tst26SteppingAction.hh"
#include "Tst26SteppingVerbose.hh"


#ifdef G4VIS_USE
#include "Tst26VisManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  //my Verbose output class
  Tst26SteppingVerbose* stepvb = new Tst26SteppingVerbose();
  G4VSteppingVerbose::SetInstance(stepvb);

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager();

  // set mandatory initialization classes
  Tst26DetectorConstruction* det = new Tst26DetectorConstruction();
  runManager->SetUserInitialization(det);
  runManager->SetUserInitialization(new Tst26PhysicsList());

#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new Tst26VisManager();
  visManager->Initialize();
#endif

  Tst26PrimaryGeneratorAction* primary = new Tst26PrimaryGeneratorAction(det);
  runManager->SetUserAction(primary);

  // set user action classes
  Tst26RunAction* RunAct = new Tst26RunAction(primary);
  runManager->SetUserAction(RunAct);
  Tst26EventAction* EvtAct = new Tst26EventAction   (RunAct);
  runManager->SetUserAction(EvtAct);
  runManager->SetUserAction(new Tst26TrackingAction(RunAct));
  runManager->SetUserAction(new Tst26SteppingAction(EvtAct));

  // get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (argc==1)   // Define UI terminal for interactive mode.
    {
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
  else           // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  delete stepvb;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

