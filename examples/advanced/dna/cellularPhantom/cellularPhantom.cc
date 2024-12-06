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
// --------------------------------------------------------------------------------
//       MONTE CARLO SIMULATION OF REALISTIC GEOMETRY FROM MICROSCOPES IMAGES
//
// Authors and contributors:
// P. Barberet, S. Incerti, N. H. Tran, L. Morelli
//
// University of Bordeaux, CNRS, LP2i, UMR5797, Gradignan, France
//
// If you use this code, please cite the following publication:
// P. Barberet et al.,
// "Monte-Carlo dosimetry on a realistic cell monolayer
// geometry exposed to alpha particles."
// Ph. Barberet et al 2012 Phys. Med. Biol. 57 2189
// doi: 110.1088/0031-9155/57/8/2189
// --------------------------------------------------------------------------------

#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4UImanager.hh"

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

int main(int argc,char** argv) {

  // Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) { ui = new G4UIExecutive(argc, argv); }

  // (Optionally) Choose the Random engine
  //G4Random::setTheEngine(new CLHEP::RanecuEngine);
  //G4Random::setTheSeed(1408);

  // Construct the default run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager();

  // Set mandatory user initialization classes
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);

  runManager->SetUserInitialization(new PhysicsList);

  // User action initialization
  runManager->SetUserInitialization(new ActionInitialization());

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  if ( ! ui ) {
    // Batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // Interactive mode
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  delete visManager;
  delete runManager;

  return 0;
}
