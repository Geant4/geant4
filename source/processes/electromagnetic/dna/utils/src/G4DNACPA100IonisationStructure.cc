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

#include "G4DNACPA100IonisationStructure.hh"

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

G4DNACPA100IonisationStructure::G4DNACPA100IonisationStructure()
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

void G4DNACPA100IonisationStructure::InitialiseGuanine()
{
  auto index = fpGuanine->GetIndex();
  ///    Guanine has 39 ionization levels

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
  fEnergyConstant[index].push_back(20.0877 * eV);
  fEnergyConstant[index].push_back(20.3321 * eV);
  fEnergyConstant[index].push_back(22.0153 * eV);
  fEnergyConstant[index].push_back(22.9048 * eV);
  fEnergyConstant[index].push_back(24.2447 * eV);
  fEnergyConstant[index].push_back(24.692 * eV);
  fEnergyConstant[index].push_back(28.2049 * eV);
  fEnergyConstant[index].push_back(32.1299 * eV);
  fEnergyConstant[index].push_back(33.2774 * eV);
  fEnergyConstant[index].push_back(33.3958 * eV);
  fEnergyConstant[index].push_back(36.6377 * eV);
  fEnergyConstant[index].push_back(37.3483 * eV);
  fEnergyConstant[index].push_back(38.3743 * eV);
  fEnergyConstant[index].push_back(305.7284 * eV);
  fEnergyConstant[index].push_back(307.4187 * eV);
  fEnergyConstant[index].push_back(307.8468 * eV);
  fEnergyConstant[index].push_back(308.9415 * eV);
  fEnergyConstant[index].push_back(309.8057 * eV);
  fEnergyConstant[index].push_back(423.1456 * eV);
  fEnergyConstant[index].push_back(423.2615 * eV);
  fEnergyConstant[index].push_back(424.5211 * eV);
  fEnergyConstant[index].push_back(425.006 * eV);
  fEnergyConstant[index].push_back(425.0315 * eV);
  fEnergyConstant[index].push_back(558.2487 * eV);
  fnLevels[index] = fEnergyConstant[index].size();
}

void G4DNACPA100IonisationStructure::InitialiseWater()
{
  auto index = fpG4_WATER->GetIndex();
  // The following values are extracted from the thesis of S. Edel,
  // Paul Sabatier University, Toulouse, France, July 7, 2006
  // Page 36
  fEnergyConstant[index].push_back(10.79 * eV);
  fEnergyConstant[index].push_back(13.39 * eV);
  fEnergyConstant[index].push_back(16.05 * eV);
  fEnergyConstant[index].push_back(32.30 * eV);
  fEnergyConstant[index].push_back(539.0 * eV);

  fUConstant[index].push_back(61.91 * eV);
  fUConstant[index].push_back(59.52 * eV);
  fUConstant[index].push_back(48.36 * eV);
  fUConstant[index].push_back(70.71 * eV);
  fUConstant[index].push_back(796.2 * eV);
  fnLevels[index] = fEnergyConstant[index].size();
}

void G4DNACPA100IonisationStructure::InitialiseDeoxyribose()
{
  auto index = fpDeoxyribose->GetIndex();
  fEnergyConstant[index].push_back(11.241 * eV);
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
  fEnergyConstant[index].push_back(18.0 * eV);
  fEnergyConstant[index].push_back(18.2696 * eV);
  fEnergyConstant[index].push_back(18.6049 * eV);
  fEnergyConstant[index].push_back(19.8378 * eV);
  fEnergyConstant[index].push_back(20.5787 * eV);
  fEnergyConstant[index].push_back(22.0146 * eV);
  fEnergyConstant[index].push_back(22.9176 * eV);
  fEnergyConstant[index].push_back(24.9005 * eV);
  fEnergyConstant[index].push_back(28.0268 * eV);
  fEnergyConstant[index].push_back(28.7417 * eV);
  fEnergyConstant[index].push_back(36.9571 * eV);
  fEnergyConstant[index].push_back(37.1881 * eV);
  fEnergyConstant[index].push_back(37.5798 * eV);
  fEnergyConstant[index].push_back(39.2622 * eV);
  fEnergyConstant[index].push_back(305.446 * eV);
  fEnergyConstant[index].push_back(306.6421 * eV);
  fEnergyConstant[index].push_back(306.8925 * eV);
  fEnergyConstant[index].push_back(307.0377 * eV);
  fEnergyConstant[index].push_back(308.5849 * eV);
  fEnergyConstant[index].push_back(559.0236 * eV);
  fEnergyConstant[index].push_back(559.3832 * eV);
  fEnergyConstant[index].push_back(559.6416 * eV);
  fEnergyConstant[index].push_back(559.7734 * eV);

  fnLevels[index] = fEnergyConstant[index].size();
}

void G4DNACPA100IonisationStructure::InitialiseCytosine()
{
  auto index = fpCytosine->GetIndex();
  ///    Cytosine has 29 ionization levels
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
  fEnergyConstant[index].push_back(20.8713 * eV);
  fEnergyConstant[index].push_back(20.9056 * eV);
  fEnergyConstant[index].push_back(24.0179 * eV);
  fEnergyConstant[index].push_back(24.5543 * eV);
  fEnergyConstant[index].push_back(29.0576 * eV);
  fEnergyConstant[index].push_back(32.0504 * eV);
  fEnergyConstant[index].push_back(34.5499 * eV);
  fEnergyConstant[index].push_back(35.5664 * eV);
  fEnergyConstant[index].push_back(38.0707 * eV);
  fEnergyConstant[index].push_back(305.7622 * eV);
  fEnergyConstant[index].push_back(307.9891 * eV);
  fEnergyConstant[index].push_back(308.674 * eV);
  fEnergyConstant[index].push_back(309.0146 * eV);
  fEnergyConstant[index].push_back(422.5331 * eV);
  fEnergyConstant[index].push_back(424.1245 * eV);
  fEnergyConstant[index].push_back(424.7781 * eV);
  fEnergyConstant[index].push_back(557.6346 * eV);
  fnLevels[index] = fEnergyConstant[index].size();
}

void G4DNACPA100IonisationStructure::InitialiseThymine()
{
  // The following values are extracted from the thesis of S. Edel,
  auto index = fpThymine->GetIndex();

  ///    THYMINE has 33 ionization levels

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
  fEnergyConstant[index].push_back(20.9944 * eV);
  fEnergyConstant[index].push_back(21.0797 * eV);
  fEnergyConstant[index].push_back(24.3676 * eV);
  fEnergyConstant[index].push_back(24.9072 * eV);
  fEnergyConstant[index].push_back(26.3981 * eV);
  fEnergyConstant[index].push_back(30.1684 * eV);
  fEnergyConstant[index].push_back(33.9007 * eV);
  fEnergyConstant[index].push_back(35.6553 * eV);
  fEnergyConstant[index].push_back(38.4935 * eV);
  fEnergyConstant[index].push_back(39.3191 * eV);
  fEnergyConstant[index].push_back(305.6808 * eV);
  fEnergyConstant[index].push_back(306.1885 * eV);
  fEnergyConstant[index].push_back(307.9374 * eV);
  fEnergyConstant[index].push_back(309.3127 * eV);
  fEnergyConstant[index].push_back(310.2121 * eV);
  fEnergyConstant[index].push_back(424.8945 * eV);
  fEnergyConstant[index].push_back(425.2178 * eV);
  fEnergyConstant[index].push_back(558.7154 * eV);
  fEnergyConstant[index].push_back(558.8106 * eV);
  fnLevels[index] = fEnergyConstant[index].size();
}

void G4DNACPA100IonisationStructure::InitialiseAdenine()
{
  auto index = fpAdenine->GetIndex();
  ///    Adenine has 35 ionization levels

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
  fEnergyConstant[index].push_back(20.2559 * eV);
  fEnergyConstant[index].push_back(21.4027 * eV);
  fEnergyConstant[index].push_back(23.0384 * eV);
  fEnergyConstant[index].push_back(23.9742 * eV);
  fEnergyConstant[index].push_back(24.479 * eV);
  fEnergyConstant[index].push_back(28.3224 * eV);
  fEnergyConstant[index].push_back(31.4796 * eV);
  fEnergyConstant[index].push_back(32.4597 * eV);
  fEnergyConstant[index].push_back(34.1263 * eV);
  fEnergyConstant[index].push_back(35.6423 * eV);
  fEnergyConstant[index].push_back(37.5026 * eV);
  fEnergyConstant[index].push_back(305.9305 * eV);
  fEnergyConstant[index].push_back(307.4532 * eV);
  fEnergyConstant[index].push_back(307.6866 * eV);
  fEnergyConstant[index].push_back(307.8303 * eV);
  fEnergyConstant[index].push_back(308.2887 * eV);
  fEnergyConstant[index].push_back(422.8443 * eV);
  fEnergyConstant[index].push_back(423.0296 * eV);
  fEnergyConstant[index].push_back(423.3954 * eV);
  fEnergyConstant[index].push_back(423.8101 * eV);
  fEnergyConstant[index].push_back(425.0749 * eV);
  fnLevels[index] = fEnergyConstant[index].size();
}

void G4DNACPA100IonisationStructure::InitialisePhosphate()
{
  auto index = fpPhosphate->GetIndex();
  ///    Phosphate has 25 ionization levels

  fEnergyConstant[index].push_back(12.9963 * eV);
  fEnergyConstant[index].push_back(12.9972 * eV);
  fEnergyConstant[index].push_back(14.3109 * eV);
  fEnergyConstant[index].push_back(15.2221 * eV);
  fEnergyConstant[index].push_back(16.0591 * eV);
  fEnergyConstant[index].push_back(16.0622 * eV);
  fEnergyConstant[index].push_back(17.6365 * eV);
  fEnergyConstant[index].push_back(17.6401 * eV);
  fEnergyConstant[index].push_back(18.8803 * eV);
  fEnergyConstant[index].push_back(20.6975 * eV);
  fEnergyConstant[index].push_back(20.7054 * eV);
  fEnergyConstant[index].push_back(24.2764 * eV);
  fEnergyConstant[index].push_back(35.6676 * eV);
  fEnergyConstant[index].push_back(38.1681 * eV);
  fEnergyConstant[index].push_back(38.1685 * eV);
  fEnergyConstant[index].push_back(40.1946 * eV);
  fEnergyConstant[index].push_back(150.138 * eV);
  fEnergyConstant[index].push_back(150.1381 * eV);
  fEnergyConstant[index].push_back(150.1414 * eV);
  fEnergyConstant[index].push_back(207.3392 * eV);
  fEnergyConstant[index].push_back(558.1119 * eV);
  fEnergyConstant[index].push_back(560.5803 * eV);
  fEnergyConstant[index].push_back(560.5808 * eV);
  fEnergyConstant[index].push_back(560.5817 * eV);
  fEnergyConstant[index].push_back(2179.592 * eV);
  fnLevels[index] = fEnergyConstant[index].size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100IonisationStructure::IonisationEnergy(const std::size_t& level,
                                                          const std::size_t& MatID)
{
  G4double ionisation = 0.;

  if (level < fnLevels[MatID]) {
    ionisation = fEnergyConstant[MatID][level];
  }
  else {
    std::ostringstream oss;
    oss << " material was not found. ";
    G4Exception("G4DNACPA100IonisationStructure::IonisationEnergy", "CPA013", FatalException,
                oss.str().c_str());
  }

  return ionisation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100IonisationStructure::UEnergy(const std::size_t& level, const std::size_t& MatID)
{
  G4double UEnergy = 0.;
  if (level < fnLevels[MatID]) {
    UEnergy = fUConstant[MatID][level];
  }
  else {
    std::ostringstream oss;
    oss << " material was not found. ";
    G4Exception("G4DNACPA100IonisationStructure::IonisationEnergy", "CPA001", FatalException,
                oss.str().c_str());
  }

  return UEnergy;
}
