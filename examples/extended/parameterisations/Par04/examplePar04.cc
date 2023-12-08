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
//
//
// --------------------------------------------------------------
//      GEANT 4 - examplePar04
// --------------------------------------------------------------
// Comments
//
// Example of a main program demonstrating inference in C++
// for fast simulation in calorimeters.
//
//-------------------------------------------------------------------
#include "FTFP_BERT.hh"                  // for FTFP_BERT
#include "Par04ActionInitialisation.hh"  // for Par04ActionInitialisation
#include "Par04DetectorConstruction.hh"  // for Par04DetectorConstruction
#include "Par04ParallelFastWorld.hh"
#include "Par04ParallelFullWorld.hh"

#include "G4EmParameters.hh"             // for G4EmParameters
#include "G4FastSimulationPhysics.hh"    // for G4FastSimulationPhysics
#include "G4HadronicProcessStore.hh"     // for G4HadronicProcessStore
#include "G4ParallelWorldPhysics.hh"
#include "G4RunManagerFactory.hh"        // for G4RunManagerFactory, G4RunMa...
#include "G4Types.hh"                    // for G4bool, G4int
#include "G4UIExecutive.hh"              // for G4UIExecutive
#include "G4UImanager.hh"                // for G4UImanager
#include "G4VisExecutive.hh"             // for G4VisExecutive

#include "G4Exception.hh"                // for G4Exception
#include "G4ExceptionSeverity.hh"        // for FatalErrorInArgument
#include "G4RunManager.hh"               // for G4RunManager
#include "G4String.hh"                   // for G4String
#include "G4VisManager.hh"               // for G4VisManager
#include "G4ios.hh"                      // for G4cout, G4endl
#include <ctime>                        // for time
#include <sstream>                       // for char_traits, operator<<, bas...
#include <string>                        // for allocator, operator+, operat...

int main(int argc, char** argv)
{
  // Macro name from arguments
  G4String batchMacroName;
  G4bool useInteractiveMode = false;
  G4int numOfThreadsOrTasks = 8;
  G4int runManagerTypeInt = 0;
  G4RunManagerType runManagerType = G4RunManagerType::Serial;
  G4String helpMsg(
    "Usage: " + G4String(argv[0]) +
    " [option(s)] \n You need to specify the mode and the macro file.\nOptions:"
    "\n\t-h\t\tdisplay this help message\n\t-m MACRO\ttriggers a batch mode "
     "executing MACRO\n\t-i\t\truns interactive mode, use it together with <-m vis*mac> macros"
    "\n\t-r\t\trun manager type (0=serial,1=MT,2=tasking)"
    "\n\t-t\t\tnumber of threads for MT mode (no change for other modes)."
    );
  if(argc < 2 ) {
    G4Exception("main", "No arguments", FatalErrorInArgument,
                ("No arguments passed to " + G4String(argv[0]) + "\n" + helpMsg)
                .c_str());
  }
  for(G4int i = 1; i < argc; ++i)
  {
    G4String argument(argv[i]);
    if(argument == "-h" || argument == "--help")
    {
      G4cout << helpMsg << G4endl;
      return 0;
    }
    else if(argument == "-m")
    {
      batchMacroName     = G4String(argv[i + 1]);
      ++i;
    }
    else if(argument == "-i")
    {
      useInteractiveMode = true;
    }
    else if(argument == "-r")
    {
      G4int tmp = atoi(argv[i + 1]);
      ++i;
      switch (tmp) {
      case 0:
        runManagerTypeInt = tmp;
        runManagerType = G4RunManagerType::Serial;
        break;
      case 1:
        runManagerTypeInt = tmp;
        runManagerType = G4RunManagerType::MTOnly;
        break;
      case 2:
        runManagerTypeInt = tmp;
        runManagerType = G4RunManagerType::Tasking;
        break;
      default:
        G4Exception("main", "Wrong Run Manager type", FatalErrorInArgument,
                "Choose 0 (serial, default), 1 (MT), 2 (tasking)");
        break;
      }
    }
    else if(argument == "-t")
    {
      numOfThreadsOrTasks = atoi(argv[i + 1]);
      ++i;
    }
    else
    {
      G4Exception("main", "Unknown argument", FatalErrorInArgument,
                  ("Unknown argument passed to " + G4String(argv[0]) + " : " +
                   argument + "\n" + helpMsg)
                    .c_str());
    }
  }

  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
  //set random seed with system time
  G4long seed = time(NULL);
  CLHEP::HepRandom::setTheSeed(seed);

  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  
  if(useInteractiveMode)
  {
    ui = new G4UIExecutive(argc, argv);
    runManagerType = G4RunManagerType::Serial;
  }

  // Initialization of default Run manager
  auto* runManager =
    G4RunManagerFactory::CreateRunManager(runManagerType);
  if(runManagerTypeInt == 1)
    runManager->SetNumberOfThreads(numOfThreadsOrTasks);
  // Detector geometry:
  auto detector = new Par04DetectorConstruction();
  auto parallelWorldFull = new Par04ParallelFullWorld("parallelWorldFullSim", detector);
  auto parallelWorldFast = new Par04ParallelFastWorld("parallelWorldFastSim", detector,
                                                      parallelWorldFull);
  detector->RegisterParallelWorld(parallelWorldFull);
  detector->RegisterParallelWorld(parallelWorldFast);
  runManager->SetUserInitialization(detector);

  // Physics list
  auto physicsList = new FTFP_BERT();
  // Add fast simulation physics
  auto fastSimulationPhysics = new G4FastSimulationPhysics();
  fastSimulationPhysics->BeVerbose();
  fastSimulationPhysics->ActivateFastSimulation("e-");
  fastSimulationPhysics->ActivateFastSimulation("e+");
  fastSimulationPhysics->ActivateFastSimulation("gamma");
  physicsList->RegisterPhysics(fastSimulationPhysics);
  // Add parallel world for readout
  physicsList->RegisterPhysics( new G4ParallelWorldPhysics("parallelWorldFullSim") );
  physicsList->RegisterPhysics( new G4ParallelWorldPhysics("parallelWorldFastSim") );
  // reduce verbosity of physics lists
  G4EmParameters::Instance()->SetVerbose(0);
  runManager->SetUserInitialization(physicsList);

  //-------------------------------
  // UserAction classes
  //-------------------------------
  runManager->SetUserInitialization(new Par04ActionInitialisation(detector, parallelWorldFull));
  //----------------
  // Visualization:
  //----------------
  G4cout << "Instantiating Visualization Manager......." << G4endl;
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if(useInteractiveMode)
  {
    if(batchMacroName.empty())
    {
      G4Exception("main", "Unknown macro name", FatalErrorInArgument,
                  ("No macro name passed to " + G4String(argv[0]))
                    .c_str());
    }
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + batchMacroName);
    ui->SessionStart();
    delete ui;
  }
  else
  {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + batchMacroName);
  }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  delete visManager;
  delete runManager;

  return 0;
}
