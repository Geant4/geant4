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
<<<<<<< HEAD
// $Id$
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
//
#include "eRositaActionInitialization.hh"
#include "eRositaDetectorConstruction.hh"
#include "eRositaPhysicsList.hh"
#include "eRositaPrimaryGeneratorAction.hh"
#include "eRositaRunAction.hh"
#include "eRositaEventAction.hh"
#include "eRositaSteppingAction.hh"
#include "eRositaSteppingVerbose.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"


int main(int argc,char** argv)
{
  // User Verbose output class
  //
  G4VSteppingVerbose* verbosity = new eRositaSteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
  
  // Run manager
  //
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  G4int nThreads = 4;
  runManager->SetNumberOfThreads(nThreads);
 
  G4cout << "***********************" << G4endl;
  G4cout << "*** Seed: " << G4Random::getTheSeed() << " ***" << G4endl;
  G4cout << "***********************" << G4endl;

  runManager->SetUserInitialization(new eRositaDetectorConstruction);
  runManager->SetUserInitialization(new eRositaPhysicsList);
  runManager->SetUserInitialization(new eRositaActionInitialization());

  // Initialize G4 kernel
  //
  runManager->Initialize();
      
  // Get the pointer to the User Interface manager
  //
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  

  if (argc != 1)   // batch mode  
  {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UImanager->ApplyCommand(command+fileName);
  }
  else // interactive mode : define visualization and UI terminal
  { 
      G4VisManager* visManager = new G4VisExecutive;
      visManager->Initialize();

      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      UImanager->ApplyCommand("/control/execute vis.mac");     
      ui->SessionStart();
      delete ui;

      delete visManager;
    }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  delete runManager;
  delete verbosity;

  return 0;
}
