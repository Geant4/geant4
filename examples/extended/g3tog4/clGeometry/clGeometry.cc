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
/// \file g3tog4/clGeometry/clGeometry.cc
/// \brief Main program of the g3tog4/clGeometry example
//
//
//
// 

// package includes
#include "G3toG4DetectorConstruction.hh"
#include "G3toG4ActionInitialization.hh"

// geant4 includes
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "FTFP_BERT.hh"
#include "G3VolTable.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr 
      << "clGeometry <call_list_file> [-m macro ] [-u UIsession] [-t nThreads]" 
      << G4endl;
    G4cerr 
      << "   note: -t option is available only for multi-threaded mode."
      << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // Evaluate arguments
  //
  G4cout << "argc " << argc << G4endl;
  if ( argc < 2 || argc > 8 ) {
    PrintUsage();
    return 1;
  }
  
  G4String inFile;
  G4String macro = "";
  G4String session;
#ifdef G4MULTITHREADED
  G4int nofThreads = 4;
#endif

  // First argument is mandatory
  inFile = argv[1];
  G4cout << "Geometry data file: " << inFile << G4endl;
  std::ifstream in(inFile);
  if ( ! in ) {
    G4cerr << "Cannot open input file \"" << inFile << "\"" << G4endl;
    return EXIT_FAILURE;
  }

  // Optional arguments
  for ( G4int i=2; i<argc; i=i+2 ) {
    G4cout << "evaluating " << argv[i] << G4endl;
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
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
    
  // Detect interactive mode (if no macro provided) and define UI session
  //
  G4UIExecutive* ui = nullptr;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  // Construct the default run manager
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(nofThreads);
#else
  G4RunManager* runManager = new G4RunManager;
#endif
    
  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetUserInitialization(new G3toG4DetectorConstruction(inFile));

  // Physics list
  runManager->SetUserInitialization(new FTFP_BERT);
    
  // User action initialization
  runManager->SetUserInitialization(new G3toG4ActionInitialization());
    
  // Initialize visualization
  //
  auto visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( macro.size() ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else  {  
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
