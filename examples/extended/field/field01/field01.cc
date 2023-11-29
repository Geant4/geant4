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

#include "F01SteppingVerbose.hh"
#include "G4RunManagerFactory.hh"

#include "F01DetectorConstruction.hh"
#include "F01ActionInitialization.hh"

#include "F01RunAction.hh"

#include "G4UImanager.hh"

// To control verbosity 
#include "G4EmParameters.hh"
#include "G4HadronicParameters.hh"

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

#include "G4TransportationParameters.hh"

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

  G4VSteppingVerbose::SetInstance(new F01SteppingVerbose);
  
  // Construct the sequential (or default) run manager
  auto* runManager =
     G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

  // G4TransportationWithMscType: fDisabled, fEnabled, fMultipleSteps
  // G4EmParameters::Instance()->SetTransportationWithMsc(G4TransportationWithMscType::fEnabled);
  
  // Set mandatory initialization classes
  //
  // Detector construction
  F01DetectorConstruction* detector = new F01DetectorConstruction();
  // detector->SetUseFSALstepper();  // Uncomment to use FSAL steppers
          
  runManager->SetUserInitialization(detector);

  // Configure the use of low thresholds for looping particles
  //  ( appropriate for typical applications using low-energy physics. )
  // auto plHelper = G4PhysicsListHelper::GetPhysicsListHelper();
  // plHelper->UseLowLooperThresholds();
  // plHelper->UseHighLooperThresholds();
  // Request a set of pre-selected values of the parameters for looping
  //    particles:
  //       - High for collider HEP applications,
  //       - Low  for 'low-E' applications, medical, ..
  // Note: If helper is used select low or high thresholds , it will overwrite
  //       values from TransportationParameters!
    
  // They are currently applied in the following order:
  // 1. Transportation Parameters - fine grained control in Transportation construction
  // 2. Physics List Helper       - impose a fixed set of new values in Transport classes
  // 3. Run Action (F01RunAction) - revise values at Start of Run
  // 4. Tracking Action           - could revise value at start of each track (not shown)
  //     Note that this also could customise by particle type, e.g. giving different values
  //     to mu-/mu+ , e-/e+ vs others)
  // If multiple are present, later methods overwrite previous ones in this list.
  
  // Physics list
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  // User action initialization
  runManager->SetUserInitialization(new F01ActionInitialization(detector));

  G4double warningE   = 10.0 * CLHEP::keV;
  G4double importantE =  0.1  * CLHEP::MeV;
  G4int    numTrials  = 30;

  G4bool useTransportParams= true; // Use the new way - Nov 2022
  
  if( useTransportParams )
  {
    auto transportParams= G4TransportationParameters::Instance();
    transportParams->SetWarningEnergy(  warningE );
    transportParams->SetImportantEnergy( importantE );
    transportParams->SetNumberOfTrials( numTrials );
    G4cout << "field01: Using G4TransportationParameters to set looper parameters."  << G4endl;
  }
  else
  {
    // Fine grained control of thresholds for looping particles
    auto runAction= new F01RunAction();
    runAction->SetWarningEnergy( warningE );
    // Looping particles with E < 10 keV will be killed after 1 step
    //   with warning.
    // Looping particles with E > 10 keV will generate a warning.
    runAction->SetImportantEnergy( importantE ); 
    runAction->SetNumberOfTrials( numTrials ); 
    // Looping particles with E > 0.1 MeV will survive for up to
    //  30 'tracking' steps, and only be killed if they still loop.

    G4cout << "field01: Using F01RunAction to set looper parameters."  << G4endl;
    runManager->SetUserAction(runAction);
  }
  
  // Note: this mechanism overwrites the thresholds established by
  //       the call to UseLowLooperThresholds() above.

  // Suppress large verbosity from EM & hadronic processes
  G4EmParameters::Instance()->SetVerbose(0);
  G4HadronicParameters::Instance()->SetVerboseLevel(0);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();

  // Initialize visualization
  //
  // G4VisExecutive can take a verbosity argument - see /vis/verbose
  G4VisManager* visManager = new G4VisExecutive("Quiet");
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
