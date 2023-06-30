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
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
//
/// \file PlasmidMolecules.cc
/// \brief Implementation of the additional Plasmid DNA molecules

#include "PlasmidMolecules.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// Deoxyribose

G4DNA_Deoxyribose* G4DNA_Deoxyribose::fDeoxyriboseInstance = 0;
G4DNA_Deoxyribose* G4DNA_Deoxyribose::Definition()
{
  
  if (fDeoxyriboseInstance != 0) return fDeoxyriboseInstance;
  const G4String name = "DNA_Deoxyribose";

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0) {
    const G4String formatedName = "DNA_Deoxy^{0}";

    G4double mass = 31.99546 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0, 0,
                                          1.7 * angstrom, 2);

    ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(0);
    ((G4MoleculeDefinition*) anInstance)->SetFormatedName(formatedName);
  }
  fDeoxyriboseInstance = static_cast<G4DNA_Deoxyribose*>(anInstance);
  return fDeoxyriboseInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// OH_Damamged_Deoxyribose

G4DNA_DamagedDeoxyriboseOH* 
  G4DNA_DamagedDeoxyriboseOH::fDamagedDeoxyriboseOHInstance = 0;
G4DNA_DamagedDeoxyriboseOH* G4DNA_DamagedDeoxyriboseOH::Definition()
{
  if (fDamagedDeoxyriboseOHInstance != 0) return fDamagedDeoxyriboseOHInstance;
  const G4String name = "DNA_DamagedDeoxyriboseOH";

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0) {
    const G4String formatedName = "DamagedDeoxyriboseOH^{0}";

    G4double mass = 31.99546 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0, 0,
                                          1.7 * angstrom, 2);

    ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(0);
    ((G4MoleculeDefinition*) anInstance)->SetFormatedName(formatedName);
  }
  fDamagedDeoxyriboseOHInstance = 
    static_cast<G4DNA_DamagedDeoxyriboseOH*>(anInstance);
  return fDamagedDeoxyriboseOHInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// H_Damamged_Deoxyribose

G4DNA_DamagedDeoxyriboseH* 
  G4DNA_DamagedDeoxyriboseH::fDamagedDeoxyriboseHInstance = 0;
G4DNA_DamagedDeoxyriboseH* G4DNA_DamagedDeoxyriboseH::Definition()
{
  if (fDamagedDeoxyriboseHInstance != 0) return fDamagedDeoxyriboseHInstance;
  const G4String name = "DNA_DamagedDeoxyriboseH";

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0) {
    const G4String formatedName = "DamagedDeoxyriboseH^{0}";

    G4double mass = 31.99546 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0, 0,
                                          1.7 * angstrom, 2);

    ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(0);
    ((G4MoleculeDefinition*) anInstance)->SetFormatedName(formatedName);
  }
  fDamagedDeoxyriboseHInstance = 
    static_cast<G4DNA_DamagedDeoxyriboseH*>(anInstance);
  return fDamagedDeoxyriboseHInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// Eaq_Damamged_Deoxyribose

G4DNA_DamagedDeoxyriboseEAQ* 
  G4DNA_DamagedDeoxyriboseEAQ::fDamagedDeoxyriboseEAQInstance = 0;
G4DNA_DamagedDeoxyriboseEAQ* G4DNA_DamagedDeoxyriboseEAQ::Definition()
{
  if (fDamagedDeoxyriboseEAQInstance != 0) return fDamagedDeoxyriboseEAQInstance;
  const G4String name = "DNA_DamagedDeoxyriboseEAQ";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0) {
    const G4String formatedName = "DamagedDeoxyriboseEAQ^{0}";

    G4double mass = 31.99546 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0, 0,
                                          1.7 * angstrom, 2);

    ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(0);
    ((G4MoleculeDefinition*) anInstance)->SetFormatedName(formatedName);
  }
  fDamagedDeoxyriboseEAQInstance = 
    static_cast<G4DNA_DamagedDeoxyriboseEAQ*>(anInstance);
  return fDamagedDeoxyriboseEAQInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
