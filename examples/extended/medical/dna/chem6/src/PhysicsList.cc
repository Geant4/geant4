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
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class

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
// Authors: W. G. Shin and S. Incerti (CENBG, France)
//

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "G4EmDNAChemistry.hh"
#include "G4EmDNAChemistry_option1.hh"
#include "G4EmDNAChemistry_option2.hh"
#include "G4EmDNAChemistry_option3.hh"
#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysics_option1.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4EmDNAPhysics_option3.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option5.hh"
#include "G4EmDNAPhysics_option6.hh"
#include "G4EmDNAPhysics_option7.hh"
#include "G4EmDNAPhysics_option8.hh"
#include "G4EmParameters.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProductionCutsTable.hh"
#include <map>
#include <functional>
#include <string>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VModularPhysicsList() {
  constexpr G4double currentDefaultCut = 0.001 * mm;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100 * eV, 1 * GeV);
  SetDefaultCutValue(currentDefaultCut);
  SetVerboseLevel(1);

  fMessenger = std::make_unique<PhysicsListMessenger>(this);

  fEmDNAPhysicsList = std::make_unique<G4EmDNAPhysics_option2>();
  fEmDNAChemistryList = std::make_unique<G4EmDNAChemistry_option3>();
  fChemDNAName = "G4EmDNAChemistry_option3";
  fPhysDNAName = "G4EmDNAPhysics_option2";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList() = default;

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

void PhysicsList::RegisterConstructor(const G4String &name) {
  // Factory maps to reduce long if-else chains and centralize constructor bindings
  using PhysFactory = std::function<std::unique_ptr<G4VPhysicsConstructor>(G4int)>;
  using ChemFactory = std::function<std::unique_ptr<G4VPhysicsConstructor>()>;

  static const std::map<std::string, PhysFactory> physFactories = {
    {"G4EmDNAPhysics", [](G4int v) { return std::make_unique<G4EmDNAPhysics>(v); }},
    {"G4EmDNAPhysics_option1", [](G4int v) { return std::make_unique<G4EmDNAPhysics_option1>(v); }},
    {"G4EmDNAPhysics_option2", [](G4int v) { return std::make_unique<G4EmDNAPhysics_option2>(v); }},
    {"G4EmDNAPhysics_option3", [](G4int v) { return std::make_unique<G4EmDNAPhysics_option3>(v); }},
    {"G4EmDNAPhysics_option4", [](G4int v) { return std::make_unique<G4EmDNAPhysics_option4>(v); }},
    {"G4EmDNAPhysics_option5", [](G4int v) { return std::make_unique<G4EmDNAPhysics_option5>(v); }},
    {"G4EmDNAPhysics_option6", [](G4int v) { return std::make_unique<G4EmDNAPhysics_option6>(v); }},
    {"G4EmDNAPhysics_option7", [](G4int v) { return std::make_unique<G4EmDNAPhysics_option7>(v); }},
    {"G4EmDNAPhysics_option8", [](G4int v) { return std::make_unique<G4EmDNAPhysics_option8>(v); }}
  };

  static const std::map<std::string, ChemFactory> chemFactories = {
    {"G4EmDNAChemistry", [] { return std::make_unique<G4EmDNAChemistry>(); }},
    {"G4EmDNAChemistry_option1", [] { return std::make_unique<G4EmDNAChemistry_option1>(); }},
    {"G4EmDNAChemistry_option2", [] { return std::make_unique<G4EmDNAChemistry_option2>(); }},
    {"G4EmDNAChemistry_option3", [] { return std::make_unique<G4EmDNAChemistry_option3>(); }}
  };

  // Try EM physics first
  if (const auto it = physFactories.find(name); it != physFactories.end()) {
    if (name == fPhysDNAName) {
      if (verboseLevel > 0) {
        G4cout << "PhysicsList: EM constructor '" << name << "' already active, no change" <<
            G4endl;
      }
      return;
    }
    if (verboseLevel > 0) {
      G4cout << "===== Register EM constructor ==== '" << name << "'" << G4endl;
    }
    fEmDNAPhysicsList = it->second(verboseLevel);
    fPhysDNAName = name;
    return;
  }

  // Then chemistry
  if (const auto it = chemFactories.find(name); it != chemFactories.end()) {
    if (name == fChemDNAName) {
      if (verboseLevel > 0) {
        G4cout << "PhysicsList: Chemistry constructor '" << name << "' already active, no change" <<
            G4endl;
      }
      return;
    }
    if (verboseLevel > 0) {
      G4cout << "===== Register Chemistry constructor ==== '" << name << "'" << G4endl;
    }
    fEmDNAChemistryList = it->second();
    if (fEmDNAChemistryList) {
      fEmDNAChemistryList->SetVerboseLevel(verboseLevel);
    }
    fChemDNAName = name;
    return;
  }

  // Unknown name
  G4cout << "PhysicsList::RegisterConstructor: <" << name << "> fails - name is not defined" <<
      G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
