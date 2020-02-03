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
/// \file genericPL.cc
/// \brief Main program of the physicslists/genericPL example
//
//
//
// -------------------------------------------------------------
//      GEANT4 Hadr00
//
//  Application demonstrating Geant4 hadronic cross sections
//
//  Author: V.Ivanchenko 20 June 2008
//
//  Modified: 
//
// -------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4GenericPhysicsList.hh"
#include "G4VModularPhysicsList.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " factory [-m macro ] [-p physListMacro ] [-u UIsession] [-t nThreads]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode." << G4endl;
    G4cerr << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) 
{
  // Evaluate arguments
  //
  if ( argc > 9 ) {
    PrintUsage();
    return 1;
  }
  
  G4String macro;
  G4String session;
  G4String physListMacro;
  G4String gdmlFileName;
#ifdef G4MULTITHREADED
  G4int nofThreads = 0;
#endif
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
    else if ( G4String(argv[i]) == "-p" ) physListMacro = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nofThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }  
  
  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  // Choose the Random engine  //choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine());

  // Construct the run manager
#ifdef G4MULTITHREADED  
  G4MTRunManager * runManager = new G4MTRunManager(); 
  if ( nofThreads > 0 ) {
    runManager->SetNumberOfThreads(nofThreads);
  }
#else
  G4RunManager * runManager = new G4RunManager(); 
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Physics List
  G4VModularPhysicsList* physList = nullptr;
  if ( physListMacro.size() ) {
    // via macro
    physList = new G4GenericPhysicsList();
    UImanager->ApplyCommand("/control/execute "+physListMacro);
  } 
  else {
    // from vector of physics cobstructor names
    std::vector<G4String>* myConstructors = new std::vector<G4String>;

    myConstructors->push_back("G4EmStandardPhysics");
    myConstructors->push_back("G4EmExtraPhysics");
    myConstructors->push_back("G4DecayPhysics");
    myConstructors->push_back("G4HadronElasticPhysics");
    myConstructors->push_back("G4HadronPhysicsFTFP_BERT");
    myConstructors->push_back("G4StoppingPhysics");
    myConstructors->push_back("G4IonPhysics");
    myConstructors->push_back("G4NeutronTrackingCut");

    physList = new G4GenericPhysicsList(myConstructors);
  }

  // Set mandatory initialization classes
  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(physList);

  // set user action classes
  runManager->SetUserInitialization(new ActionInitialization("genericPL"));

  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  if ( macro.size() ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else {
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac"); 
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
