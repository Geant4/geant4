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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file main.cc
/// \brief Chem3 example

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4DNAChemistryManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
 * WARNING : Geant4 was initially not intended for this kind of application
 * This code is delivered as a prototype
 * We will be happy to hear from you, do not hesitate to send your feedback and
 * communicate on the difficulties you may encounter
 * The user interface may change in the next releases since a reiteration of the
 * code has started
 */

int main(int argc, char** argv)
{
  // Construct the default run manager

  G4RunManager * runManager;
  G4bool isMT(true);

  if (isMT == true)
  { // use it but just with one thread (test purpose)
#ifdef G4MULTITHREADED
  runManager= new G4MTRunManager;
  ((G4MTRunManager*)runManager)->SetNumberOfThreads(2);
#else
    runManager = new G4RunManager();
#endif

  }
  else
  {
    runManager = new G4RunManager;
  }

  // activate chemistery
  G4bool chemistryFlag(true);
  G4DNAChemistryManager::Instance()->SetChemistryActivation(chemistryFlag);

  // Set mandatory user initialization classes
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserInitialization(detector);

  G4DNAChemistryManager::Instance()->InitializeMaster();
  runManager->SetUserInitialization(new ActionInitialization());

  // Initialize G4 kernel
  runManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  if (argc != 1)
  {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }
  else
  {
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute vis.mac");
#endif
    ui->SessionStart();
    delete ui;
#else
    UImanager->ApplyCommand("/control/execute beam.in");
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
