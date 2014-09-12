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
// $Id$
//
/// \file testAnalysis.cc
/// \brief Main program of the  testAnalysis

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "ApplicationParameters.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT_HP.hh"

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

using namespace ApplicationParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " testAnalysis "
           << "    [-m macro ]" << G4endl
           << "    [-all    on|off]" << G4endl
           << "    [-hn     on|off]" << G4endl
           << "    [-pn     on|off]" << G4endl
           << "    [-h1     on|off]" << G4endl
           << "    [-h2     on|off]" << G4endl
           << "    [-h3     on|off]" << G4endl
           << "    [-p1     on|off]" << G4endl
           << "    [-p2     on|off]" << G4endl
           << "    [-ntuple on|off]" << G4endl
           << "    [-read   on|off]" << G4endl
           << "    [-write  on|off]" << G4endl
           << G4endl;
  }

  G4bool GetOption(G4String option) {
    if      ( option == "on" )  return true;
    else if ( option == "off" ) return false;
    else  {
      PrintUsage();
      return 1;
    }
  }

  void TestMacros() {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String h1Macro = "h1.mac";
    G4String h2Macro = "h2.mac";
    G4String h3Macro = "h3.mac";
    G4String p1Macro = "p1.mac";
    G4String p2Macro = "p2.mac";
    if ( TestH1 ) UImanager->ApplyCommand(command + h1Macro);
    if ( TestH2 ) UImanager->ApplyCommand(command + h2Macro);
    if ( TestH3 ) UImanager->ApplyCommand(command + h3Macro);
    if ( TestP1 ) UImanager->ApplyCommand(command + p1Macro);
    if ( TestP2 ) UImanager->ApplyCommand(command + p2Macro);
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  G4String macro;
  for ( G4int i=1; i<argc; i=i+2 ) {
    //G4cout << "... testing " << G4String(argv[i]) << G4endl;
    if      ( G4String(argv[i]) == "-m" )   macro = argv[i+1];
    else if ( G4String(argv[i]) == "-all" ) SetTestAll(GetOption(argv[i+1]));
    else if ( G4String(argv[i]) == "-hn" )  SetTestHn(GetOption(argv[i+1]));
    else if ( G4String(argv[i]) == "-pn" )  SetTestPn(GetOption(argv[i+1]));
    else if ( G4String(argv[i]) == "-h1" )  TestH1 = GetOption(argv[i+1]);
    else if ( G4String(argv[i]) == "-h2" )  TestH2 = GetOption(argv[i+1]);
    else if ( G4String(argv[i]) == "-h3" )  TestH3 = GetOption(argv[i+1]);
    else if ( G4String(argv[i]) == "-p1" )  TestP1 = GetOption(argv[i+1]);
    else if ( G4String(argv[i]) == "-p2" )  TestP2 = GetOption(argv[i+1]);
    else if ( G4String(argv[i]) == "-ntuple" )  TestNtuple = GetOption(argv[i+1]);
    else if ( G4String(argv[i]) == "-read" )    TestRead = GetOption(argv[i+1]);
    else if ( G4String(argv[i]) == "-write" )   TestWrite = GetOption(argv[i+1]);
    else {
      PrintUsage();
      return 1;
    }
  }  
  
  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  //
  DetectorConstruction* detConstruction = new DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  G4VModularPhysicsList* physicsList = new QGSP_BERT_HP;
  runManager->SetUserInitialization(physicsList);
    
  ActionInitialization* actionInitialization
     = new ActionInitialization();
  runManager->SetUserInitialization(actionInitialization);

  // Initialize G4 kernel
  //
  runManager->Initialize();
  
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( macro.size() ) {
    // batch mode
    TestMacros();
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else  {  
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac"); 
#else
    UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
    if (ui->IsGUI())
      UImanager->ApplyCommand("/control/execute gui.mac");
    TestMacros();
    ui->SessionStart();
    delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
