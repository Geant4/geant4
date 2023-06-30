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
#include "PhysicsList.hh"
#include "DetectorConstruction.hh"
#include "EmDNAChemistry.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4PhysicsConstructorRegistry.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList(DetectorConstruction *pDetector)
    : G4VModularPhysicsList(),
      fEmDNAPhysicsList(new G4EmDNAPhysics_option2(verboseLevel)),
      fEmDNAChemistryList(new EmDNAChemistry) {
  G4double currentDefaultCut = 0.01 * mm;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100 * eV,
                                                                  1 * GeV);
  SetDefaultCutValue(currentDefaultCut);
  fEmDNAChemistryList->SetChemistryWorld(*(pDetector->GetChemistryWorld()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle() {
  if (fEmDNAPhysicsList != nullptr) {
    fEmDNAPhysicsList->ConstructParticle();
  }
  if (fEmDNAChemistryList != nullptr) {
    fEmDNAChemistryList->ConstructParticle();
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess() {
  AddTransportation();
  if (fEmDNAPhysicsList != nullptr) {
    fEmDNAPhysicsList->ConstructProcess();
  }
  if (fEmDNAChemistryList != nullptr) {
    fEmDNAChemistryList->ConstructProcess();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......