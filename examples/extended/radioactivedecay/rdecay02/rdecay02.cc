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
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "exrdmDetectorConstruction.hh"
#include "exrdmPhysicsList.hh"
#include "exrdmEventAction.hh"
#include "exrdmRunAction.hh"
#include "exrdmSteppingAction.hh"
#include "exrdmPrimaryGeneratorAction.hh"
#include "exrdmAnalysisManager.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv)
{
  // random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // Creation of the analysis manager
  exrdmAnalysisManager::getInstance();

  // set mandatory initialization classes

  exrdmDetectorConstruction* Detector = new exrdmDetectorConstruction;
  runManager->SetUserInitialization(Detector);
  runManager->SetUserInitialization(new exrdmPhysicsList);

  // set mandatory user action class
  runManager->SetUserAction(new exrdmPrimaryGeneratorAction);
  runManager->SetUserAction(new exrdmRunAction);
  runManager->SetUserAction(new exrdmEventAction);
  runManager->SetUserAction(new exrdmSteppingAction);

#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // Initialize G4 kernel
  // do this at run time so the geometry/physics can be changed
  //  runManager->Initialize();

  // get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if (argc == 1)   // Define UI session for interactive mode.
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
  exrdmAnalysisManager::dispose();
  delete runManager;

  return 0;
}








