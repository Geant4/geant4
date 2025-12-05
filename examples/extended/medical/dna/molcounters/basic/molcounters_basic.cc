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
/// \file molcounters_basic.cc
/// \brief Main program of the molcounters/basic example

// The `molcounters` example(s) are provided as part of Geant4-DNA
// and any report or published result obtained using it shall cite
// the respective Geant4-DNA collaboration publications.
//
// Reports or results obtained using the spatially-aware `MoleculeCounter`
// provided in this example, shall further cite:
//
// Velten & Tomé, Radiation Physics and Chemistry, 2023 (10.1016/j.radphyschem.2023.111194)
//
//
// Author: Christian Velten (2025)
//

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include "G4DNAChemistryManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

int main(int argc, char** argv)
{
  G4UIExecutive* ui = 0;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  auto* runManager = G4RunManagerFactory::CreateRunManager();

  // Set mandatory initialization classes
  runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserInitialization(new DetectorConstruction);
  runManager->SetUserInitialization(new ActionInitialization);

  // get the pointer to the User Interface manager
  auto UI = G4UImanager::GetUIpointer();

  if (argc > 1)  // batch mode
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command + fileName);
  }
  else  // define visualization and UI terminal for interactive mode
  {
    UI->ApplyCommand("/control/execute simple_sbs.in");
    delete ui;
  }

  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
