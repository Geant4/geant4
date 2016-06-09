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
// ----------------------------------------------------------------------------
//                 GEANT4 - GAMMAKNIFE example
// ----------------------------------------------------------------------------
// AUTHORS:
// G. Cuttone (a), J. Pipek (b) F.Romano* (a,c), M.G.Sabini (d)
// 
// PAST AUTHORS:
// G.A.P. Cirrone (a), G.Russo (e), M.Russo (a)
//
// (a) Laboratori Nazionali del Sud - INFN, Catania, Italy
// (b) Faculty of Nuclear Sciences and Physical Engineering, Czech Technical University, Czech Republic
// (c) Centro Studi e Ricerche e Museo Storico della Fisica E.Fermi, Roma, Italy
// (d) Dipartimento di Immagini, Ospedale Cannizzaro, Catania, Italy
// (e) Fondazione Istituto San Raffaele G.Giglio, Cefalù (Palermo), Italy
//
//
// *Corresponding author, email to francesco.romano@lns.infn.it
// ----------------------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif
#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif
#include "GammaKnifeDetectorConstruction.hh"
#include "GammaKnifePhysicsList.hh"
#include "GammaKnifePrimaryGeneratorAction.hh"
#include "GammaKnifeRunAction.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4UImessenger.hh"
#include "G4ScoringManager.hh"
#include "globals.hh"
#include "GammaKnifeController.hh"
#include "G4PhysListFactory.hh"

#include <ctime>

int main(int argc ,char ** argv)
{
  
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  G4int seconds =  time(NULL);
  CLHEP::HepRandom::setTheSeed(seconds);
  
  G4RunManager* runManager = new G4RunManager;
  G4ScoringManager::GetScoringManager(); // This enables scoring

  // Initialize the geometry
  GammaKnifeDetectorConstruction* detector = new GammaKnifeDetectorConstruction();
  runManager -> SetUserInitialization(detector);
 
	
  // Initialize the physics 
  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = 0;
  G4String physName = "";

  // Physics List name defined via environment variable
  char* path = getenv("PHYSLIST");
  if (path) { physName = G4String(path); }

  if(physName != "" && factory.IsReferencePhysList(physName))
    {
      phys = factory.GetReferencePhysList(physName);
    } 

  if(!phys) { phys = new GammaKnifePhysicsList(); }

  runManager->SetUserInitialization(phys);

  // Initialize the primary particles  
  runManager -> SetUserAction(new GammaKnifePrimaryGeneratorAction());

  GammaKnifeController* controller = new GammaKnifeController( detector );
  controller->ReadFile("MachineAngle.in"); // pre-load default

  // Optional UserActions: run, event, stepping
  GammaKnifeRunAction* pRunAction = new GammaKnifeRunAction();
  runManager -> SetUserAction(pRunAction);

  // Initialize G4 kernel
  //
  runManager->Initialize();
  
#ifdef G4VIS_USE
  // Visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager -> Initialize();
#endif

  // Get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  
 
  if (argc!=1) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute defaultMacro.mac"); 
#else
    UImanager->ApplyCommand("/control/execute batch.mac"); 
#endif
    ui->SessionStart();
    delete ui;
#endif
  }

  // Job termination
#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;
  delete controller;

  return 0;
}
