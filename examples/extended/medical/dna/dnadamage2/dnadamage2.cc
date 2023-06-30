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
// chem6 example is derived from chem4 and chem5 examples
//
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// J. Appl. Phys. 125 (2019) 104301
// Med. Phys. 45 (2018) e722-e739
// J. Comput. Phys. 274 (2014) 841-882
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157-178
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//          Jose Ramos-Mendez and Bruce Faddegon (UCSF, US)
//
/// \file dnadamage2.cc
/// \brief DnaDamage2 example

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#include "G4RunManagerFactory.hh"

#include "G4DNAChemistryManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

/*
 * WARNING : Geant4 was initially not intended for this kind of application
 * This code is delivered as a prototype
 * We will be happy to hear from you, do not hesitate to send your feedback
 * and communicate on the difficulties you may encounter
 * The user interface may change in the next releases since a reiteration of
 * the code has started
 */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4int fSeed = 1234;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

int main(int argc, char** argv)
{
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  auto* runManager= G4RunManagerFactory::CreateRunManager();

  // Set mandatory initialization classes
  DetectorConstruction* fpDetector  = new DetectorConstruction();
  ActionInitialization* fpActionIni = new ActionInitialization();

  runManager->SetUserInitialization(new PhysicsList());
  runManager->SetUserInitialization(fpDetector);
  runManager->SetUserInitialization(fpActionIni);

  //get the pointer to the User Interface manager
  G4UImanager* UI  = G4UImanager::GetUIpointer();
  G4VisManager* vM = new G4VisExecutive;

  G4String fileName = "";
  G4String command = "/control/execute ";

  if (argc == 1) {
    vM->Initialize();
    G4Random::setTheSeed(fSeed);
    UI->ApplyCommand("/control/execute init_vis.in");
    ui->SessionStart();
    delete ui;
  }

  else if (argc == 2) // batch mode
  {
    fileName = argv[1];
  }

  else if (argc > 2) 
  {
    fileName = argv[1];
    fSeed = atoi(argv[2]);
  }

  if (argc > 1) {
    G4Random::setTheSeed(fSeed);
    UI->ApplyCommand(command+fileName);
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
  delete runManager;
  delete vM;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
