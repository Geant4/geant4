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
/// \file exoticphysics/monopole/monopole.cc
/// \brief Main program of the exoticphysics/monopole example
//
<<<<<<< HEAD
// $Id: monopole.cc 92500 2015-09-02 07:26:32Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include "G4MonopolePhysics.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "RunAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr
      << " Usage: " << G4endl
      << " monopole [-m macro ] [-s setupMonopole] [-t nThreads]" << G4endl
      << "   Note: " << G4endl
      << "    -s should be followed by a composed string, eg. \'1 0 100 GeV\'" << G4endl
      << "    -t option is for multi-threaded mode." << G4endl
      << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 7 ) {
    PrintUsage();
    return 1;
  }

  G4String macro;
  G4String setupMonopole;
  G4int nThreads = 1;
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-s" ) setupMonopole = argv[i+1];
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
    else {
      PrintUsage();
      return 1;
    }
  }

  // Construct the default run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  if ( nThreads > 0 ) {
    runManager->SetNumberOfThreads(nThreads);
  }
  G4cout << "===== Example is started with "
         <<  runManager->GetNumberOfThreads() << " threads =====" << G4endl;

  // Instantiate G4UIExecutive if interactive
  G4UIExecutive* ui = nullptr;
  if ( macro.empty() ) {
    ui = new G4UIExecutive(argc, argv);
  }

  //get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  //create physicsList
  // Physics List is defined via environment variable PHYSLIST
  G4PhysListFactory factory;
<<<<<<< HEAD
  G4VModularPhysicsList* phys = factory.ReferencePhysList(); 
 
=======
  G4VModularPhysicsList* phys = factory.GetReferencePhysList("FTFP_BERT");

>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  // monopole physics is added
  G4MonopolePhysics * theMonopole = new G4MonopolePhysics();

  // Setup monopole
  if ( setupMonopole.size() )  {
    UImanager->ApplyCommand("/control/verbose 1");
    UImanager->ApplyCommand("/monopole/setup " + setupMonopole);
  }

  // regsiter monopole physics
  phys->RegisterPhysics(theMonopole);

  runManager->SetUserInitialization(phys);

  // visualization manager
  G4VisManager* visManager = nullptr;

  // set detector construction
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  DetectorConstruction* det = new DetectorConstruction();

  runManager->SetUserInitialization(det);  
  runManager->SetUserInitialization(phys);

<<<<<<< HEAD
  PrimaryGeneratorAction* kin = new PrimaryGeneratorAction(det);
  runManager->SetUserAction(kin);

  //user action classes
  RunAction* run;

  runManager->SetUserAction(run = new RunAction(det, kin));
  runManager->SetUserAction(new TrackingAction(run));
  runManager->SetUserAction(new SteppingAction(run));

  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);    
    }
  else
    {  // interactive mode : define UI session
#ifdef G4VIS_USE
      //visualization manager
      visManager = new G4VisExecutive;
      visManager->Initialize();
#endif
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      ui->SessionStart();
      delete ui;
#endif
    }
=======
  // Process macro or start UI session
  //
  if ( macro.size() ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else  {
    // interactive mode : define UI session
    visManager = new G4VisExecutive();
    visManager->Initialize();
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

  //job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
 
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
