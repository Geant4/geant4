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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Models come from
// M. Bug et al, Rad. Phys and Chem. 130, 459-479 (2017)

#include "G4DNAPTBIonisationStructure.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

G4DNAPTBIonisationStructure::G4DNAPTBIonisationStructure()
{
  fpTHF = G4Material::GetMaterial("THF", false);
  fpPY = G4Material::GetMaterial("PY", false);
  fpPU = G4Material::GetMaterial("PU", false);
  fpTMP = G4Material::GetMaterial("TMP", false);
  fpG4_WATER = G4Material::GetMaterial("G4_WATER", false);
  fpBackbone_THF = G4Material::GetMaterial("backbone_THF", false);
  fpCytosine_PY = G4Material::GetMaterial("cytosine_PY", false);
  fpThymine_PY = G4Material::GetMaterial("thymine_PY", false);
  fpAdenine_PU = G4Material::GetMaterial("adenine_PU", false);
  fpBackbone_TMP = G4Material::GetMaterial("backbone_TMP", false);
  fpGuanine_PU = G4Material::GetMaterial("guanine_PU", false);
  fpN2 = G4Material::GetMaterial("N2", false);

  // MPietrzak
  if (fpN2 != nullptr) {
    auto index = fpN2->GetIndex();
    energyConstant[index].push_back(15.58 * eV);
    energyConstant[index].push_back(17.07 * eV);
    energyConstant[index].push_back(21.00 * eV);
    energyConstant[index].push_back(41.72 * eV);
  }

  // MPietrzak
  if (fpG4_WATER != nullptr) {
    auto index = fpG4_WATER->GetIndex();
    energyConstant[index].push_back(10.79 * eV);
    energyConstant[index].push_back(13.39 * eV);
    energyConstant[index].push_back(16.05 * eV);
    energyConstant[index].push_back(32.30 * eV);
    energyConstant[index].push_back(539.0 * eV);
  }
  if (fpTHF != nullptr) {
    auto index = fpTHF->GetIndex();
    energyConstant[index].push_back(9.74 * eV);
    energyConstant[index].push_back(12.31 * eV);
    energyConstant[index].push_back(12.99 * eV);
    energyConstant[index].push_back(13.57 * eV);
    energyConstant[index].push_back(13.60 * eV);
    energyConstant[index].push_back(15.11 * eV);
    energyConstant[index].push_back(15.97 * eV);
    energyConstant[index].push_back(16.28 * eV);
    energyConstant[index].push_back(18.19 * eV);
    energyConstant[index].push_back(18.69 * eV);
    energyConstant[index].push_back(22.14 * eV);
    energyConstant[index].push_back(22.25 * eV);
    energyConstant[index].push_back(27.21 * eV);
    energyConstant[index].push_back(28.97 * eV);
    energyConstant[index].push_back(36.97 * eV);
    energyConstant[index].push_back(305.07 * eV);
    energyConstant[index].push_back(305.08 * eV);
    energyConstant[index].push_back(306.17 * eV);
    energyConstant[index].push_back(306.17 * eV);
    energyConstant[index].push_back(557.94 * eV);
  }

  if (fpPY != nullptr) {
    auto index = fpPY->GetIndex();
    energyConstant[index].push_back(9.73 * eV);
    energyConstant[index].push_back(10.96 * eV);
    energyConstant[index].push_back(11.54 * eV);
    energyConstant[index].push_back(12.58 * eV);
    energyConstant[index].push_back(15.96 * eV);
    energyConstant[index].push_back(16.27 * eV);
    energyConstant[index].push_back(16.53 * eV);
    energyConstant[index].push_back(17.98 * eV);
    energyConstant[index].push_back(19.37 * eV);
    energyConstant[index].push_back(20.52 * eV);
    energyConstant[index].push_back(24.55 * eV);
    energyConstant[index].push_back(24.64 * eV);
    energyConstant[index].push_back(29.75 * eV);
    energyConstant[index].push_back(33.02 * eV);
    energyConstant[index].push_back(36.57 * eV);
    energyConstant[index].push_back(305.92 * eV);
    energyConstant[index].push_back(307.09 * eV);
    energyConstant[index].push_back(307.09 * eV);
    energyConstant[index].push_back(307.52 * eV);
    energyConstant[index].push_back(423.44 * eV);
    energyConstant[index].push_back(423.44 * eV);
  }

  if (fpPU != nullptr) {
    auto index = fpPU->GetIndex();
    energyConstant[index].push_back(9.58 * eV);
    energyConstant[index].push_back(10.57 * eV);
    energyConstant[index].push_back(10.97 * eV);
    energyConstant[index].push_back(12.22 * eV);
    energyConstant[index].push_back(12.92 * eV);
    energyConstant[index].push_back(13.44 * eV);
    energyConstant[index].push_back(15.05 * eV);
    energyConstant[index].push_back(16.56 * eV);
    energyConstant[index].push_back(17.18 * eV);
    energyConstant[index].push_back(17.88 * eV);
    energyConstant[index].push_back(17.90 * eV);
    energyConstant[index].push_back(19.11 * eV);
    energyConstant[index].push_back(20.09 * eV);
    energyConstant[index].push_back(21.70 * eV);
    energyConstant[index].push_back(23.52 * eV);
    energyConstant[index].push_back(24.35 * eV);
    energyConstant[index].push_back(25.41 * eV);
    energyConstant[index].push_back(29.34 * eV);
    energyConstant[index].push_back(32.44 * eV);
    energyConstant[index].push_back(33.67 * eV);
    energyConstant[index].push_back(36.26 * eV);
    energyConstant[index].push_back(38.22 * eV);
    energyConstant[index].push_back(306.53 * eV);
    energyConstant[index].push_back(307.19 * eV);
    energyConstant[index].push_back(307.64 * eV);
    energyConstant[index].push_back(308.14 * eV);
    energyConstant[index].push_back(308.17 * eV);
    energyConstant[index].push_back(423.31 * eV);
    energyConstant[index].push_back(423.43 * eV);
    energyConstant[index].push_back(423.64 * eV);
    energyConstant[index].push_back(423.98 * eV);
  }

  if (fpTMP != nullptr) {
    auto index = fpTMP->GetIndex();
    energyConstant[index].push_back(10.81 * eV);
    energyConstant[index].push_back(10.81 * eV);
    energyConstant[index].push_back(12.90 * eV);
    energyConstant[index].push_back(13.32 * eV);
    energyConstant[index].push_back(13.32 * eV);
    energyConstant[index].push_back(13.59 * eV);
    energyConstant[index].push_back(14.33 * eV);
    energyConstant[index].push_back(14.33 * eV);
    energyConstant[index].push_back(15.90 * eV);
    energyConstant[index].push_back(17.09 * eV);
    energyConstant[index].push_back(17.09 * eV);
    energyConstant[index].push_back(17.13 * eV);
    energyConstant[index].push_back(17.85 * eV);
    energyConstant[index].push_back(17.85 * eV);
    energyConstant[index].push_back(18.44 * eV);
    energyConstant[index].push_back(19.37 * eV);
    energyConstant[index].push_back(19.37 * eV);
    energyConstant[index].push_back(21.40 * eV);
    energyConstant[index].push_back(26.20 * eV);
    energyConstant[index].push_back(26.20 * eV);
    energyConstant[index].push_back(27.43 * eV);
    energyConstant[index].push_back(35.23 * eV);
    energyConstant[index].push_back(37.67 * eV);
    energyConstant[index].push_back(37.67 * eV);
    energyConstant[index].push_back(39.64 * eV);
    energyConstant[index].push_back(152.42 * eV);
    energyConstant[index].push_back(152.42 * eV);
    energyConstant[index].push_back(152.44 * eV);
    energyConstant[index].push_back(209.59 * eV);
    energyConstant[index].push_back(306.92 * eV);
    energyConstant[index].push_back(306.92 * eV);
    energyConstant[index].push_back(306.92 * eV);
    energyConstant[index].push_back(557.34 * eV);
    energyConstant[index].push_back(559.40 * eV);
    energyConstant[index].push_back(559.40 * eV);
    energyConstant[index].push_back(559.41 * eV);
    energyConstant[index].push_back(2178.05 * eV);
  }

  for (const auto& [index, levels] : energyConstant) {
    nLevels[index] = (G4int)levels.size();
  }
}

G4double G4DNAPTBIonisationStructure::IonisationEnergy(G4int level, const size_t& materialID)
{
  size_t matNameModif = ReplaceMaterial(materialID);

  // check if the material exist in the map
  if (energyConstant.find(matNameModif) == energyConstant.end()) {
    std::ostringstream oss;
    oss << "Material name was not found in energyConstantMap. Problematic material is: "
        << materialID;
    G4Exception(
      "G4DNAPTBIonisationStructure::IonisationEnergy", "em0002", FatalException, oss.str().c_str());
  }

  G4double ionisation = 0.;

  if (level >= 0 && level < nLevels[matNameModif]) ionisation = energyConstant[matNameModif][level];

  return ionisation;
}

G4int G4DNAPTBIonisationStructure::NumberOfLevels(const size_t& materialID)
{
  auto matNameModif = ReplaceMaterial(materialID);

  // check if the material exist in the map
  if (nLevels.find(matNameModif) == nLevels.end()) {
    std::ostringstream oss;
    oss << "Material name was not found in energyConstantMap. Problematic material is: "
        << matNameModif;
    G4Exception(
      "G4DNAPTBIonisationStructure::NumberOfLevels", "em0002", FatalException, oss.str().c_str());
  }

  return nLevels[matNameModif];
}

size_t G4DNAPTBIonisationStructure::ReplaceMaterial(const size_t& materialID)
{
  if (fpBackbone_THF != nullptr && materialID == fpBackbone_THF->GetIndex()) {
    return fpTHF->GetIndex();
  }
  else if (fpBackbone_TMP != nullptr && materialID == fpBackbone_TMP->GetIndex()) {
    return fpTMP->GetIndex();
  }
  else if (fpAdenine_PU != nullptr && materialID == fpAdenine_PU->GetIndex()) {
    return fpPU->GetIndex();
  }
  else if (fpGuanine_PU != nullptr && materialID == fpGuanine_PU->GetIndex()) {
    return fpPU->GetIndex();
  }
  else if (fpThymine_PY != nullptr && materialID == fpThymine_PY->GetIndex()) {
    return fpPY->GetIndex();
  }
  else if (fpCytosine_PY != nullptr && materialID == fpCytosine_PY->GetIndex()) {
    return fpPY->GetIndex();
  }
  return materialID;
}
