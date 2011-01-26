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
// $Id: test46.cc,v 1.3 2010-12-21 19:43:08 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------
//      GEANT4 test46
//
//  Simplified combined CMS calorimeter 
//
//  Authors: A.Ivanchenko and V.Ivanchenko
//           September 2008
//
//  Modified: 
//
// -------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "StackingAction.hh"
#include "G4Timer.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  G4Timer timer;
  timer.Start();

  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());

  //Construct the default run manager
  G4RunManager * runManager = new G4RunManager();

  //set mandatory initialization classes
  runManager->SetUserInitialization(new DetectorConstruction());

  // Physics List name defined via environment variable
  // or the default test46 name is used
  G4String physName = "";
  char* path = getenv("PHYSLIST");
  if(path) { physName = G4String(path); }
  if("" == physName) { physName = "FTFP_BERT_EMV"; }

  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = factory.GetReferencePhysList(physName);

  runManager->SetUserInitialization(phys);

  runManager->SetUserAction(new PrimaryGeneratorAction());

  //set user action classes
  runManager->SetUserAction(new RunAction());
  runManager->SetUserAction(new EventAction());
  runManager->SetUserAction(new StackingAction());

  // get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (argc==1)   // Define UI terminal for interactive mode
    {
#ifdef G4VIS_USE
      //visualization manager
      G4VisExecutive* visManager = new G4VisExecutive;
      visManager->Initialize();
#endif

#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
      ui->SessionStart();
      delete ui;
#endif

#ifdef G4VIS_USE
     delete visManager;
#endif
    }
  else           // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  //job termination
  delete runManager;
  timer.Stop();
  G4cout << timer << std::endl;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
