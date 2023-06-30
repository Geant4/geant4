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
/// \file ScavengerMolecules.cc
/// \brief Implementation of the Background Scavenber chemical species

#include "ScavengerMolecules.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4DMSO* G4DMSO::fDMSOInstance = 0;

G4DMSO* G4DMSO::Definition()
{
  if (fDMSOInstance != 0) return fDMSOInstance;
  const G4String name = "DMSO^0";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0) {
    const G4String formatedName = "DMSO^{0}";

    G4double mass = 31.99546 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 2.4e-9 * (m * m / s), 0, 0,
                                          1.7 * angstrom, 2);

    ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(0); // not implemented
    ((G4MoleculeDefinition*) anInstance)->SetFormatedName(formatedName);
  }
  fDMSOInstance = static_cast<G4DMSO*>(anInstance);
  return fDMSOInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4OxygenB* G4OxygenB::fOxygenBInstance = 0;

G4OxygenB* G4OxygenB::Definition()
{
  if (fOxygenBInstance != 0) return fOxygenBInstance;
  const G4String name = "Oxygen(B)";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0) {
    const G4String formatedName = "Oxygen(B)^{0}";

    G4double mass = 31.99546 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 2.4e-9 * (m * m / s), 0, 0,
                                          1.7 * angstrom, 2);

    ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(0); // not implemented
    ((G4MoleculeDefinition*) anInstance)->SetFormatedName(formatedName);
  }
  fOxygenBInstance = static_cast<G4OxygenB*>(anInstance);
  return fOxygenBInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
