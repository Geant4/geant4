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
// Based on the work of M. Terrissol and M. C. Bordage
//
// Users are requested to cite the following papers:
// - M. Terrissol, A. Baudre, Radiat. Prot. Dosim. 31 (1990) 175-177
// - M.C. Bordage, J. Bordes, S. Edel, M. Terrissol, X. Franceries,
//   M. Bardies, N. Lampe, S. Incerti, Phys. Med. 32 (2016) 1833-1840
//
// Authors of this class:
// M.C. Bordage, M. Terrissol, S. Edel, J. Bordes, S. Incerti
//
// 15.01.2014: creation
//
//
// Modified for Adenine material by S. Zein on 20.04.2021
// Based on the study by S. Zein et. al. Nucl. Inst. Meth. B 488 (2021) 70-82

#include "G4DNACPA100ExcitationStructure.hh"

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
G4DNACPA100ExcitationStructure::G4DNACPA100ExcitationStructure()
{
  fpGuanine = G4Material::GetMaterial("G4_GUANINE", false);
  fpG4_WATER = G4Material::GetMaterial("G4_WATER", false);
  fpDeoxyribose = G4Material::GetMaterial("G4_DEOXYRIBOSE", false);
  fpCytosine = G4Material::GetMaterial("G4_CYTOSINE", false);
  fpThymine = G4Material::GetMaterial("G4_THYMINE", false);
  fpAdenine = G4Material::GetMaterial("G4_ADENINE", false);
  fpPhosphate = G4Material::GetMaterial("G4_PHOSPHORIC_ACID", false);

  if (fpGuanine != nullptr) {
    InitialiseGuanine();
  }

  if (fpG4_WATER != nullptr) {
    InitialiseWater();
  }

  if (fpDeoxyribose != nullptr) {
    InitialiseDeoxyribose();
  }

  if (fpCytosine != nullptr) {
    InitialiseCytosine();
  }

  if (fpThymine != nullptr) {
    InitialiseThymine();
  }

  if (fpAdenine != nullptr) {
    InitialiseAdenine();
  }

  if (fpPhosphate != nullptr) {
    InitialisePhosphate();
  }
}

void G4DNACPA100ExcitationStructure::InitialiseGuanine()
{
  auto index = fpGuanine->GetIndex();
  fEnergyConstant[index].push_back(8.2315 * eV);
  fEnergyConstant[index].push_back(11.0928 * eV);
  fEnergyConstant[index].push_back(11.5984 * eV);
  fEnergyConstant[index].push_back(11.7906 * eV);
  fEnergyConstant[index].push_back(11.9382 * eV);
  fEnergyConstant[index].push_back(12.4424 * eV);
  fEnergyConstant[index].push_back(13.3581 * eV);
  fEnergyConstant[index].push_back(15.1381 * eV);
  fEnergyConstant[index].push_back(16.4059 * eV);
  fEnergyConstant[index].push_back(16.496 * eV);
  fEnergyConstant[index].push_back(16.8457 * eV);
  fEnergyConstant[index].push_back(17.297 * eV);
  fEnergyConstant[index].push_back(18.0608 * eV);
  fEnergyConstant[index].push_back(18.441 * eV);
  fEnergyConstant[index].push_back(19.2414 * eV);

  fnLevels[index] = fEnergyConstant[index].size();
}
void G4DNACPA100ExcitationStructure::InitialiseWater()
{
  auto index = fpG4_WATER->GetIndex();
  // The following values are extracted from the thesis of S. Edel,
  // Paul Sabatier University, Toulouse, France, July 7, 2006
  // Page 36
  fEnergyConstant[index].push_back(8.17 * eV);
  fEnergyConstant[index].push_back(10.13 * eV);
  fEnergyConstant[index].push_back(11.31 * eV);
  fEnergyConstant[index].push_back(12.91 * eV);
  fEnergyConstant[index].push_back(14.50 * eV);

  fUConstant[index].push_back(61.91 * eV);
  fUConstant[index].push_back(59.52 * eV);
  fUConstant[index].push_back(48.36 * eV);
  fUConstant[index].push_back(70.71 * eV);
  fUConstant[index].push_back(796.2 * eV);

  fnLevels[index] = fEnergyConstant[index].size();
}
void G4DNACPA100ExcitationStructure::InitialiseDeoxyribose()
{
  auto index = fpDeoxyribose->GetIndex();
  fEnergyConstant[index].push_back(11.2410 * eV);
  fEnergyConstant[index].push_back(11.7927 * eV);
  fEnergyConstant[index].push_back(12.6579 * eV);
  fEnergyConstant[index].push_back(12.8163 * eV);
  fEnergyConstant[index].push_back(13.3238 * eV);
  fEnergyConstant[index].push_back(13.9487 * eV);
  fEnergyConstant[index].push_back(14.4374 * eV);
  fEnergyConstant[index].push_back(14.7433 * eV);
  fEnergyConstant[index].push_back(15.0818 * eV);
  fEnergyConstant[index].push_back(15.6112 * eV);
  fEnergyConstant[index].push_back(16.0547 * eV);
  fEnergyConstant[index].push_back(16.8319 * eV);
  fEnergyConstant[index].push_back(17.4294 * eV);
  fEnergyConstant[index].push_back(18.0000 * eV);
  fEnergyConstant[index].push_back(18.2696 * eV);
  fEnergyConstant[index].push_back(18.6049 * eV);
  fEnergyConstant[index].push_back(19.8378 * eV);

  fnLevels[index] = fEnergyConstant[index].size();
}
void G4DNACPA100ExcitationStructure::InitialiseCytosine()
{
  auto index = fpCytosine->GetIndex();
  // The following values are extracted from the thesis of S. Edel,
  fEnergyConstant[index].push_back(9.3222 * eV);
  fEnergyConstant[index].push_back(10.4601 * eV);
  fEnergyConstant[index].push_back(11.3044 * eV);
  fEnergyConstant[index].push_back(11.9986 * eV);
  fEnergyConstant[index].push_back(13.4528 * eV);
  fEnergyConstant[index].push_back(14.7371 * eV);
  fEnergyConstant[index].push_back(16.2286 * eV);
  fEnergyConstant[index].push_back(16.5877 * eV);
  fEnergyConstant[index].push_back(17.0741 * eV);
  fEnergyConstant[index].push_back(17.1875 * eV);
  fEnergyConstant[index].push_back(18.638 * eV);
  fEnergyConstant[index].push_back(19.6884 * eV);
  fnLevels[index] = fEnergyConstant[index].size();
}
void G4DNACPA100ExcitationStructure::InitialiseThymine()
{
  // The following values are extracted from the thesis of S. Edel,
  auto index = fpThymine->GetIndex();

  fEnergyConstant[index].push_back(9.639 * eV);
  fEnergyConstant[index].push_back(11.8278 * eV);
  fEnergyConstant[index].push_back(12.0876 * eV);
  fEnergyConstant[index].push_back(12.9656 * eV);
  fEnergyConstant[index].push_back(13.9555 * eV);
  fEnergyConstant[index].push_back(15.0774 * eV);
  fEnergyConstant[index].push_back(15.4078 * eV);
  fEnergyConstant[index].push_back(15.4689 * eV);
  fEnergyConstant[index].push_back(16.1964 * eV);
  fEnergyConstant[index].push_back(16.8955 * eV);
  fEnergyConstant[index].push_back(17.5018 * eV);
  fEnergyConstant[index].push_back(18.2979 * eV);
  fEnergyConstant[index].push_back(18.4495 * eV);
  fEnergyConstant[index].push_back(19.3186 * eV);
  fnLevels[index] = fEnergyConstant[index].size();
}
void G4DNACPA100ExcitationStructure::InitialiseAdenine()
{
  auto index = fpAdenine->GetIndex();
  // The following values are extracted from the thesis of S. Edel,

  fEnergyConstant[index].push_back(8.5114 * eV);
  fEnergyConstant[index].push_back(10.1294 * eV);
  fEnergyConstant[index].push_back(11.0606 * eV);
  fEnergyConstant[index].push_back(11.5849 * eV);
  fEnergyConstant[index].push_back(12.1533 * eV);
  fEnergyConstant[index].push_back(13.356 * eV);
  fEnergyConstant[index].push_back(13.6554 * eV);
  fEnergyConstant[index].push_back(15.3296 * eV);
  fEnergyConstant[index].push_back(16.179 * eV);
  fEnergyConstant[index].push_back(16.7676 * eV);
  fEnergyConstant[index].push_back(17.3489 * eV);
  fEnergyConstant[index].push_back(17.5568 * eV);
  fEnergyConstant[index].push_back(18.554 * eV);
  fEnergyConstant[index].push_back(19.0866 * eV);
  fnLevels[index] = fEnergyConstant[index].size();
}
void G4DNACPA100ExcitationStructure::InitialisePhosphate()
{
  // The following values are extracted from the thesis of S. Edel,
  auto index = fpPhosphate->GetIndex();
  fEnergyConstant[index].push_back(12.9963 * eV);
  fEnergyConstant[index].push_back(12.9972 * eV);
  fEnergyConstant[index].push_back(14.3109 * eV);
  fEnergyConstant[index].push_back(15.2221 * eV);
  fEnergyConstant[index].push_back(16.0591 * eV);
  fEnergyConstant[index].push_back(16.0622 * eV);
  fEnergyConstant[index].push_back(17.6365 * eV);
  fEnergyConstant[index].push_back(17.6401 * eV);
  fEnergyConstant[index].push_back(18.8803 * eV);
  fnLevels[index] = fEnergyConstant[index].size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ExcitationStructure::ExcitationEnergy(const std::size_t& level,
                                                          const std::size_t& MatID)
{
  G4double ionisation = 0.;

  if (level < fnLevels[MatID]) {
    ionisation = fEnergyConstant[MatID][level];
  }
  else {
    std::ostringstream oss;
    oss << " material was not found. ";
    G4Exception("G4DNACPA100ExcitationStructure::ExcitationEnergy", "CPA001", FatalException,
                oss.str().c_str());
  }

  return ionisation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ExcitationStructure::UEnergy(const std::size_t& level, const std::size_t& MatID)
{
  G4double UEnergy = 0.;
  if (level < fnLevels[MatID]) {
    UEnergy = fUConstant[MatID][level];
  }
  else {
    std::ostringstream oss;
    oss << " material was not found. ";
    G4Exception("G4DNACPA100ExcitationStructure::ExcitationEnergy", "CPA001", FatalException,
                oss.str().c_str());
  }

  return UEnergy;
}
