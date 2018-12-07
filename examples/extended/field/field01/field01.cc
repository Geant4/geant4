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
/// \file field/field01/field01.cc
/// \brief Main program of the field/field01 example
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"

#ifdef G4MULTITHREADED
// #define USE_MULTITHREADED
#endif

#ifdef   USE_MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "F01SteppingVerbose.hh"
#include "G4RunManager.hh"
#endif

#include "F01DetectorConstruction.hh"
#include "F01ActionInitialization.hh"

#include "F01RunAction.hh"

#include "G4UImanager.hh"

#include "G4EmParameters.hh"
#include "G4HadronicProcessStore.hh"

#include "G4PhysicsListHelper.hh"

#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"
#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

// For Printing statistic from Transporation process(es)
#include "G4Electron.hh"
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Instantiate G4UIExecutive if there are no arguments (interactive mode)
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //
#ifdef USE_MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
#else
  G4VSteppingVerbose::SetInstance(new F01SteppingVerbose);
  G4RunManager * runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  //
  // Detector construction
  F01DetectorConstruction* detector = new F01DetectorConstruction();
  // detector->SetUseFSALstepper();  // Uncomment to use FSAL steppers
          
  runManager->SetUserInitialization(detector);

  // Configure the use of low thresholds for looping particles
  //  ( appropriate for typical applications using low-energy physics. )
  auto plHelper = G4PhysicsListHelper::GetPhysicsListHelper();
  plHelper->UseLowLooperThresholds();
  // Request a set of pre-selected values of the parameters for looping
  //  particles
  
  // Physics list
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  // User action initialization
  runManager->SetUserInitialization(new F01ActionInitialization(detector));

  // Fine grained control of thresholds for looping particles
  auto runAction= new F01RunAction();
  runAction->SetWarningEnergy(   10.0 * CLHEP::keV );
              // Looping particles with E < 10 keV will be killed after 1 step
              //   with warning.
              // Looping particles with E > 10 keV will generate a warning.
  runAction->SetImportantEnergy( 0.1  * CLHEP::MeV );
  runAction->SetNumberOfTrials( 30 );
              // Looping particles with E > 0.1 MeV will survive for up to
              //  30 'tracking' steps, and only be killed if they still loop.
  // Note: this mechanism overwrites the thresholds established by
  //       the call to UseLowLooperThresholds() above.
  
  runManager->SetUserAction(runAction);

  // Suppress large verbosity from EM & hadronic processes
  G4EmParameters::Instance()->SetVerbose(-1);
  G4HadronicProcessStore::Instance()->SetVerbose(0);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (!ui)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  else
    {  // interactive mode : define UI session
     UImanager->ApplyCommand("/control/execute init_vis.mac");
     if (ui->IsGUI())
       UImanager->ApplyCommand("/control/execute gui.mac");
     ui->SessionStart();
     delete ui;
    }

  // Statistics of tracks killed by G4Transportation are currently
  //  printed in the RunAction's EndOfEvent action.
  // ( Eventually a summary could be provided here instead or as well. )

  delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
