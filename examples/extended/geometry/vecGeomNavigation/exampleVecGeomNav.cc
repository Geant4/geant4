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
//  Author: J. Apostolakis, S. Wenzel,  2018-2021
//
//  Started from FullCMS Geant4 application by Mihaly Novak  (2017)
//---------------------------------------------------------------------

#include <iostream>
#include <iomanip>

#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4UIterminal.hh"

#include "G4VModularPhysicsList.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"

#include "Randomize.hh"

#include "VG01DetectorConstruction.hh"
#include "VG01ActionInitialization.hh"
#include "VG01SteppingVerboseWithDir.hh"

// For interactivity
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

static G4bool       parUseVecGeom = true;
static G4bool       parCompareG4  = false;
static G4bool       parInteractive = false;
static std::string  parMacroFileName = "";
static std::string  parGDMLFile = "TestNTST.gdml";

void GetInputArguments(int argc, char** argv);
void PrintUsage();

int main(int argc, char** argv) {
  //
  // get input arguments
  GetInputArguments(argc, argv);
  G4cout<< " ========== Running exampleVecGeomNav ================ " << G4endl
        << "   GDML geometry file          =  " << parGDMLFile       << G4endl
        << "   Geant4 macro                =  " << parMacroFileName  << G4endl
        << "   Use VecGeom (VG) navigation =  " << parUseVecGeom     << G4endl
        << "   Compare G4 vs VG navigation =  " << parCompareG4      << G4endl
        << " ===================================================== " << G4endl;

  // Use custom stepping verbosity
  G4VSteppingVerbose::SetInstance( new VG01SteppingVerboseWithDir() );       

  // Construct the run manager
  //
  G4RunManager* runManager 
      = G4RunManagerFactory::CreateRunManager( G4RunManagerType::Serial);
  //  or G4RunManagerType::Default to get Task or Multithreading

  // set mandatory initialization classes
  //
  // 1. Detector construction
  //
  VG01DetectorConstruction* detector = new VG01DetectorConstruction;
  detector->SetGDMLFileName(parGDMLFile);
  detector->SetUseVecGeom(parUseVecGeom);
  runManager->SetUserInitialization(detector);

  // 2. Physics list
  //
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  // 3. User action
  //  
  runManager->SetUserInitialization(new VG01ActionInitialization());

  // 4. Run the simulation in batch mode, except if macroFile == "-"
  //  
  G4UImanager* UImgr = G4UImanager::GetUIpointer();
  G4String command = "/control/execute ";
  if( parMacroFileName != "-")
     UImgr->ApplyCommand(command+parMacroFileName);

  // 5. Run the simulation in Interactive mode if requested ( flag: -i )
  //
  if( parInteractive )
  {
     G4UIExecutive* uiExec = 0;
     uiExec = new G4UIExecutive(argc, argv);

     // Initialize visualization
     //
     G4VisManager* visManager = new G4VisExecutive;
     // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
     // G4VisManager* visManager = new G4VisExecutive("Quiet");
     visManager->Initialize();

     // UImgr->ApplyCommand("/control/execute init_vis.mac");

     // interactive mode
     uiExec->SessionStart();

     // Cleanup
     delete uiExec;
     delete visManager;
  } 
  else 
  {

    // Print out the final random number - for batch only runs
    G4cout << G4endl
           << " ================================================================= " << G4endl
           << " Final random number = " << G4UniformRand() << G4endl
           << " ================================================================= " << G4endl
           << G4endl;
  }
  //

  // Delete the RunManager
  delete runManager;
  return 0;
}

void horizontal_line(char c)
{
  std::cout <<"\n " << std::setw(100) << std::setfill(c) 
            << "" << std::setfill(' ') << std::endl;
}

void PrintUsage() {
  horizontal_line('=');
  std::cout << "  Geant4 application to demonstrate interface to VecGeom Navigation.    \n"
            << std::endl
            << "  Two modes: \n\n"
            << "   * 1 parameter this is treated as Geant4 macro file \n"
            << " \n"
            << "   * Multiple Parameters: \n"
            << "      at least one of the following: \n"
            << "       -m :   the standard Geant4 macro file \n"
            << "       -i :   interactive (after batch, if any) \n"
            << "      optionally one of the following: \n"
            << "       -v :   flag  ==> run using VecGeom navigation (default). \n"
            << "       -o :   flag  ==> run using Geant4  navigation. \n"
            << "       -c :   flag  ==> compare VecGeom and Geant4 navigation"             
            << " (and report differences.) \n"
            << "      and other(s): \n"
            << "       -g :   GDML file with geometry \n"
            << "\n"
            << std::endl;
  horizontal_line('=');
}

void GetInputArguments(int argc, char** argv) {
  // process arguments
  if (argc == 1) {
    PrintUsage();
    exit(0);
  }
  if (argc == 2) {
     parMacroFileName = argv[1];
     G4cout << " argc = 2  -- Filename = " << parMacroFileName << G4endl;
     return;
  }

  // Adapted from examples/basic/B4/exampleB4a.cc
  for ( G4int i=1; i<argc; ++i ) {
    if  ( G4String(argv[i]) == "-m" ) {
      if( ++i < argc ) {
        parMacroFileName = argv[i];
        G4cout << " arg-parsing:  Macro file name= " << parMacroFileName << G4endl;
      }
      else{
         G4cerr << " Parse Error: '-m' cannot be last argument.  Use it : -m <filename>" << G4endl;
         PrintUsage();
         exit(1);
      }
    }
    else if ( G4String(argv[i]) == "-g" ) {
      if( ++i < argc ) { parGDMLFile      = argv[i];
         G4cout << " arg-parsing:  GDML file name= " << parGDMLFile << G4endl;
      }
      else{
         G4cerr << " Parse Error: '-m' cannot be last argument.  Use it : -m <filename>" << G4endl;
         PrintUsage();
         exit(1);
      }
    }
    else if ( G4String(argv[i]) == "-v" ) { parUseVecGeom = true; }
    else if ( G4String(argv[i]) == "-o" ) { parUseVecGeom = false; }
    else if ( G4String(argv[i]) == "-c" ) { parCompareG4 = true;  parUseVecGeom= true; }
    else {
      G4cerr << "  Unknown argument : " << argv[i] << G4endl;
      PrintUsage();
      exit(1);
    }
  }

  // check if mandatory Geant4 macro file was provided
  if (parMacroFileName=="" && !parInteractive ) {
     G4cout << "  *** ERROR : either interactive mode or a Geant4 macro file is required. " << G4endl;
     PrintUsage();
     exit(-1);
  }
}
