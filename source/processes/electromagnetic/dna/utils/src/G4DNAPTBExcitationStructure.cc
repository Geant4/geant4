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

#include "G4DNAPTBExcitationStructure.hh"

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

G4DNAPTBExcitationStructure::G4DNAPTBExcitationStructure()
{
  fpN2 = G4Material::GetMaterial("N2", false);

  //  taken directly from PTra code by MPietrzak
  if (fpN2 != nullptr) {
    auto index = fpN2->GetIndex();
    energyConstant[index].push_back(1.85 * eV);
    energyConstant[index].push_back(2.15 * eV);
    energyConstant[index].push_back(8.00 * eV);
    energyConstant[index].push_back(8.50 * eV);
    energyConstant[index].push_back(8.60 * eV);
    energyConstant[index].push_back(11.05 * eV);
    energyConstant[index].push_back(11.79 * eV);
    energyConstant[index].push_back(11.90 * eV);
    energyConstant[index].push_back(12.25 * eV);
    energyConstant[index].push_back(12.50 * eV);
    energyConstant[index].push_back(13.01 * eV);
    energyConstant[index].push_back(13.19 * eV);
    energyConstant[index].push_back(13.30 * eV);
    energyConstant[index].push_back(14.33 * eV);
    energyConstant[index].push_back(14.84 * eV);
    energyConstant[index].push_back(15.18 * eV);
    energyConstant[index].push_back(15.70 * eV);
    energyConstant[index].push_back(15.75 * eV);
    energyConstant[index].push_back(15.86 * eV);
    energyConstant[index].push_back(17.36 * eV);
    energyConstant[index].push_back(17.95 * eV);
    energyConstant[index].push_back(19.77 * eV);
    energyConstant[index].push_back(20.79 * eV);
    energyConstant[index].push_back(20.87 * eV);
    energyConstant[index].push_back(22.27 * eV);
    energyConstant[index].push_back(22.83 * eV);
    energyConstant[index].push_back(37.19 * eV);
    energyConstant[index].push_back(38.67 * eV);
    energyConstant[index].push_back(39.23 * eV);
  }

  for (const auto& [index, levels] : energyConstant) {
    nExcLevels[index] = (G4int)levels.size();
  }
}

G4double G4DNAPTBExcitationStructure::ExcitationEnergy(
  const G4int& ExcLevel, const size_t& materialID)
{
  size_t matNameModif = ReplaceMaterial(materialID);

  // check if the material exist in the map
  if (energyConstant.find(matNameModif) == energyConstant.end()) {
    std::ostringstream oss;
    oss << "Material name was not found in energyConstantMap. Problematic material is: "
        << matNameModif;
    G4Exception(
      "G4DNAPTBExcitationStructure::ExcitationEnergy", "em0002", FatalException, oss.str().c_str());
  }

  G4double excitation = 0.;

  if (ExcLevel >= 0 && ExcLevel < nExcLevels[matNameModif])
    excitation = energyConstant[matNameModif][ExcLevel];

  return excitation;
}

G4int G4DNAPTBExcitationStructure::NumberOfExcLevels(const size_t& matID)
{
  auto matNameModif = ReplaceMaterial(matID);

  // check if the material exist in the map
  if (nExcLevels.find(matNameModif) == nExcLevels.end()) {
    std::ostringstream oss;
    oss << "Material name was not found in energyConstantMap. Problematic material is: "
        << matNameModif;
    G4Exception("G4DNAPTBNDExcitationStructure::NumberOfExcLevels", "em0002", FatalException,
      oss.str().c_str());
  }

  return nExcLevels[matNameModif];
}

size_t G4DNAPTBExcitationStructure::ReplaceMaterial(const size_t& materialID)
{
  auto output = materialID;
  auto G4_N2 = G4Material::GetMaterial("G4_N2", false)->GetIndex();
  if (materialID == G4_N2) {
    output = fpN2->GetIndex();
  }

  return output;
}
