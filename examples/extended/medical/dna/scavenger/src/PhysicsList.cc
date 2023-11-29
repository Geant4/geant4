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
/// \file scavenger/src/PhysicsList.cc
/// \brief Implementation of the scavenger::PhysicsList class

#include "PhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "EmDNAChemistry.hh"
#include "G4EmDNAPhysics_option2.hh"

// // Change the physicsList
// #include "G4EmDNAPhysics.hh"
// #include "G4EmDNAPhysics_option1.hh"
// #include "G4EmDNAPhysics_option3.hh"
// #include "G4EmDNAPhysics_option4.hh"
// #include "G4EmDNAPhysics_option5.hh"
// #include "G4EmDNAPhysics_option6.hh"
// #include "G4EmDNAPhysics_option7.hh"
// #include "G4EmDNAPhysics_option8.hh"

#include "G4PhysicsConstructorRegistry.hh"
#include "G4EmParameters.hh"

namespace scavenger
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
  : G4VModularPhysicsList(),
    fpEmDNAPhysicsList(new G4EmDNAPhysics_option2(verboseLevel)),
    fpEmDNAChemistryList(new EmDNAChemistry()) {
  G4double currentDefaultCut = 1. * nanometer;
  G4ProductionCutsTable::GetProductionCutsTable()->
    SetEnergyRange(100 * eV, 1 * GeV);
  SetDefaultCutValue(currentDefaultCut);
  SetVerboseLevel(1);
  fpEmDNAPhysicsList->SetVerboseLevel(verboseLevel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle() {
  if (fpEmDNAPhysicsList) {
    fpEmDNAPhysicsList->ConstructParticle();
  }
  if (fpEmDNAChemistryList) {
    fpEmDNAChemistryList->ConstructParticle();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess() {
  AddTransportation();
  if (fpEmDNAPhysicsList) {
    fpEmDNAPhysicsList->ConstructProcess();
  }
  if (fpEmDNAChemistryList) {
    fpEmDNAChemistryList->ConstructProcess();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}