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
// CPA100 elastic model class for electrons
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
// Based on the study by S. Zein et. al. Nucl. Inst. Meth. B 488 (2021) 70-82
// 1/2/2023 : Hoang added modification

#include "G4DNACPA100ElasticModel.hh"

#include "G4DNAMaterialManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4SystemOfUnits.hh"
using namespace std;
G4DNACPA100ElasticModel::G4DNACPA100ElasticModel(const G4ParticleDefinition*, const G4String& nam)
  : G4VDNAModel(nam, "all")
{
  fpGuanine = G4Material::GetMaterial("G4_GUANINE", false);
  fpG4_WATER = G4Material::GetMaterial("G4_WATER", false);
  fpDeoxyribose = G4Material::GetMaterial("G4_DEOXYRIBOSE", false);
  fpCytosine = G4Material::GetMaterial("G4_CYTOSINE", false);
  fpThymine = G4Material::GetMaterial("G4_THYMINE", false);
  fpAdenine = G4Material::GetMaterial("G4_ADENINE", false);
  fpPhosphate = G4Material::GetMaterial("G4_PHOSPHORIC_ACID", false);
  fpParticle = G4Electron::ElectronDefinition();
}

void G4DNACPA100ElasticModel::Initialise(const G4ParticleDefinition* p,
                                         const G4DataVector& /*cuts*/)
{
  if (isInitialised) {
    return;
  }
  if (verboseLevel > 3) {
    G4cout << "Calling G4DNACPA100ExcitationModel::Initialise()" << G4endl;
  }

  if (!G4DNAMaterialManager::Instance()->IsLocked()) {
    if (p != fpParticle) {
      std::ostringstream oss;
      oss << " Model is not applied for this particle " << p->GetParticleName();
      G4Exception("G4DNACPA100ElasticModel::G4DNACPA100ElasticModel", "CPA001", FatalException,
                  oss.str().c_str());
    }

    char* path = getenv("G4LEDATA");

    if (!path) {
      G4Exception("G4DNACPA100ElasticModel::Initialise", "em0006", FatalException,
                  "G4LEDATA environment variable not set.");
      return;
    }

    std::size_t index;
    if (fpG4_WATER != nullptr) {
      index = fpG4_WATER->GetIndex();
      fLevels[index] = 1.214e-4;
      AddCrossSectionData(index, p, "dna/sigma_elastic_e_cpa100",
                          "dna/sigmadiff_cumulated_elastic_e_cpa100", 1e-20 * m * m);
      SetLowELimit(index, p, 11. * eV);
      SetHighELimit(index, p, 255955. * eV);
    }
    if (fpGuanine != nullptr) {
      index = fpGuanine->GetIndex();
      fLevels[index] = 1.4504480e-05;
      AddCrossSectionData(index, p, "dna/sigma_elastic_e_cpa100_guanine",
                          "dna/sigmadiff_cumulated_elastic_e_cpa100_guanine", 1 * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpDeoxyribose != nullptr) {
      index = fpDeoxyribose->GetIndex();
      fLevels[index] = 1.6343100e-05;
      AddCrossSectionData(index, p, "dna/sigma_elastic_e_cpa100_deoxyribose",
                          "dna/sigmadiff_cumulated_elastic_e_cpa100_deoxyribose", 1 * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpCytosine != nullptr) {
      index = fpCytosine->GetIndex();
      fLevels[index] = 1.9729660e-05;
      AddCrossSectionData(index, p, "dna/sigma_elastic_e_cpa100_cytosine",
                          "dna/sigmadiff_cumulated_elastic_e_cpa100_cytosine", 1 * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpThymine != nullptr) {
      index = fpThymine->GetIndex();
      fLevels[index] = 1.7381300e-05;
      AddCrossSectionData(index, p, "dna/sigma_elastic_e_cpa100_thymine",
                          "dna/sigmadiff_cumulated_elastic_e_cpa100_thymine", 1 * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpAdenine != nullptr) {
      index = fpAdenine->GetIndex();
      fLevels[index] = 1.6221800e-05;
      AddCrossSectionData(index, p, "dna/sigma_elastic_e_cpa100_adenine",
                          "dna/sigmadiff_cumulated_elastic_e_cpa100_adenine", 1 * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpPhosphate != nullptr) {
      index = fpPhosphate->GetIndex();
      fLevels[index] = 2.2369600e-05;
      AddCrossSectionData(index, p, "dna/sigma_elastic_e_cpa100_phosphoric_acid",
                          "dna/sigmadiff_cumulated_elastic_e_cpa100_phosphoric_acid", 1 * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }

    // Load data
    LoadCrossSectionData(p);
    G4DNAMaterialManager::Instance()->SetMasterDataModel(DNAModelType::fDNAElastics, this);
    fpModelData = this;
  }
  else {
    auto dataModel = dynamic_cast<G4DNACPA100ElasticModel*>(
      G4DNAMaterialManager::Instance()->GetModel(DNAModelType::fDNAElastics));
    if (dataModel == nullptr) {
      G4cout << "G4DNACPA100ElasticModel::CrossSectionPerVolume:: not good modelData" << G4endl;
      G4Exception("G4DNACPA100ElasticModel::CrossSectionPerVolume", "em004", FatalException,
                  "no modelData is registered");
    }
    else {
      fpModelData = dataModel;
    }
  }

  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ElasticModel::CrossSectionPerVolume(const G4Material* pMaterial,
                                                        const G4ParticleDefinition* p,
                                                        G4double ekin, G4double, G4double)
{
  // Get the name of the current particle
  const G4String& particleName = p->GetParticleName();
  auto materialID = pMaterial->GetIndex();

  // set killBelowEnergy value for current material
  fKillBelowEnergy = fpModelData->GetLowELimit(materialID, p);

  G4double sigma = 0.;

  if (ekin < fpModelData->GetHighELimit(materialID, p)) {
    if (ekin < fKillBelowEnergy) {
      return DBL_MAX;
    }

    auto tableData = fpModelData->GetData();

    if ((*tableData)[materialID][p] == nullptr) {
      G4Exception("G4DNACPA100ElasticModel::CrossSectionPerVolume", "em00236", FatalException,
                  "No model is registered");
    }
    sigma = (*tableData)[materialID][p]->FindValue(ekin);
  }

  if (verboseLevel > 2) {
    auto MolDensity =
      (*G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(pMaterial))[materialID];
    G4cout << "__________________________________" << G4endl;
    G4cout << "°°° G4DNACPA100ElasticModel - XS INFO START" << G4endl;
    G4cout << "°°° Kinetic energy(eV)=" << ekin / eV << " particle : " << particleName << G4endl;
    G4cout << "°°° lowLim (eV) = " << GetLowELimit(materialID, p) / eV
           << " highLim (eV) : " << GetHighELimit(materialID, p) / eV << G4endl;
    G4cout << "°°° Materials = " << (*G4Material::GetMaterialTable())[materialID]->GetName()
           << G4endl;
    G4cout << "°°° Cross section per molecule (cm^2)=" << sigma / cm / cm << G4endl;
    G4cout << "°°° Cross section per Phosphate molecule (cm^-1)=" << sigma * MolDensity / (1. / cm)
           << G4endl;
    G4cout << "°°° G4DNACPA100ElasticModel - XS INFO END" << G4endl;
  }

  auto MolDensity =
    (*G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(pMaterial))[materialID];
  return sigma * MolDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNACPA100ElasticModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                                                const G4MaterialCutsCouple* couple,
                                                const G4DynamicParticle* aDynamicElectron, G4double,
                                                G4double)
{
  G4double electronEnergy0 = aDynamicElectron->GetKineticEnergy();
  auto materialID = couple->GetMaterial()->GetIndex();
  auto p = aDynamicElectron->GetParticleDefinition();

  if (p != fpParticle) {
    G4Exception("G4DNACPA100ElasticModel::SampleSecondaries", "em00436", FatalException,
                "This particle is not applied for this model");
  }
  if (electronEnergy0 < fKillBelowEnergy) {
    return;
  }
  G4double cosTheta = fpModelData->RandomizeCosTheta(electronEnergy0, materialID);
  G4double phi = 2. * CLHEP::pi * G4UniformRand();

  const G4ThreeVector& zVers = aDynamicElectron->GetMomentumDirection();

  G4double CT1, ST1, CF1, SF1, CT2, ST2, CF2, SF2;
  G4double sinTheta = std::sqrt(1 - cosTheta * cosTheta);

  CT1 = zVers.z();
  ST1 = std::sqrt(1. - CT1 * CT1);

  if (ST1 != 0)
    CF1 = zVers.x() / ST1;
  else
    CF1 = std::cos(2. * CLHEP::pi * G4UniformRand());
  if (ST1 != 0)
    SF1 = zVers.y() / ST1;
  else
    SF1 = std::sqrt(1. - CF1 * CF1);

  G4double A3, A4, A5, A2, A1;

  A3 = sinTheta * std::cos(phi);
  A4 = A3 * CT1 + ST1 * cosTheta;
  A5 = sinTheta * std::sin(phi);
  A2 = A4 * SF1 + A5 * CF1;
  A1 = A4 * CF1 - A5 * SF1;

  CT2 = CT1 * cosTheta - ST1 * A3;
  ST2 = std::sqrt(1. - CT2 * CT2);

  if (ST2 == 0) ST2 = 1E-6;
  CF2 = A1 / ST2;
  SF2 = A2 / ST2;
  G4ThreeVector zPrimeVers(ST2 * CF2, ST2 * SF2, CT2);

  fParticleChangeForGamma->ProposeMomentumDirection(zPrimeVers.unit());

  auto EnergyDeposit = fpModelData->GetElasticLevel(materialID) * (1. - cosTheta) * electronEnergy0;
  fParticleChangeForGamma->ProposeLocalEnergyDeposit(EnergyDeposit);
  if (statCode) {
    fParticleChangeForGamma->SetProposedKineticEnergy(electronEnergy0);
  }
  else {
    auto newEnergy = electronEnergy0 - EnergyDeposit;
    fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
  }
}

G4double G4DNACPA100ElasticModel::Theta(const G4ParticleDefinition* p, G4double k,
                                        G4double integrDiff, const std::size_t& materialID)
{
  G4double theta, valueT1, valueT2, valueE21, valueE22, valueE12, valueE11;
  G4double xs11 = 0;
  G4double xs12 = 0;
  G4double xs21 = 0;
  G4double xs22 = 0;
  if (p == G4Electron::ElectronDefinition()) {
    if (k == tValuesVec[materialID][p].back()) {
      k = k * (1. - 1e-12);
    }
    auto t2 =
      std::upper_bound(tValuesVec[materialID][p].begin(), tValuesVec[materialID][p].end(), k);
    auto t1 = t2 - 1;

    auto e12 = std::upper_bound(eValuesVect[materialID][p][(*t1)].begin(),
                                eValuesVect[materialID][p][(*t1)].end(), integrDiff);
    auto e11 = e12 - 1;

    auto e22 = std::upper_bound(eValuesVect[materialID][p][(*t2)].begin(),
                                eValuesVect[materialID][p][(*t2)].end(), integrDiff);
    auto e21 = e22 - 1;

    valueT1 = *t1;
    valueT2 = *t2;
    valueE21 = *e21;
    valueE22 = *e22;
    valueE12 = *e12;
    valueE11 = *e11;

    xs11 = diffCrossSectionData[materialID][p][valueT1][valueE11];
    xs12 = diffCrossSectionData[materialID][p][valueT1][valueE12];
    xs21 = diffCrossSectionData[materialID][p][valueT2][valueE21];
    xs22 = diffCrossSectionData[materialID][p][valueT2][valueE22];
  }

  if (xs11 == 0 && xs12 == 0 && xs21 == 0 && xs22 == 0) {
    return (0.);
  }

  theta = QuadInterpolator(valueE11, valueE12, valueE21, valueE22, xs11, xs12, xs21, xs22, valueT1,
                           valueT2, k, integrDiff);

  return theta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ElasticModel::LinLogInterpolate(G4double e1, G4double e2, G4double e,
                                                    G4double xs1, G4double xs2)
{
  G4double d1 = std::log(xs1);
  G4double d2 = std::log(xs2);
  G4double value = std::exp(d1 + (d2 - d1) * (e - e1) / (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ElasticModel::LinLinInterpolate(G4double e1, G4double e2, G4double e,
                                                    G4double xs1, G4double xs2)
{
  G4double d1 = xs1;
  G4double d2 = xs2;
  G4double value = (d1 + (d2 - d1) * (e - e1) / (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ElasticModel::LogLogInterpolate(G4double e1, G4double e2, G4double e,
                                                    G4double xs1, G4double xs2)
{
  G4double a = (std::log10(xs2) - std::log10(xs1)) / (std::log10(e2) - std::log10(e1));
  G4double b = std::log10(xs2) - a * std::log10(e2);
  G4double sigma = a * std::log10(e) + b;
  G4double value = (std::pow(10., sigma));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ElasticModel::QuadInterpolator(G4double e11, G4double e12, G4double e21,
                                                   G4double e22, G4double xs11, G4double xs12,
                                                   G4double xs21, G4double xs22, G4double t1,
                                                   G4double t2, G4double t, G4double e)
{
  // Log-Log
  /*
    G4double interpolatedvalue1 = LogLogInterpolate(e11, e12, e, xs11, xs12);
    G4double interpolatedvalue2 = LogLogInterpolate(e21, e22, e, xs21, xs22);
    G4double value = LogLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);


    // Lin-Log
    G4double interpolatedvalue1 = LinLogInterpolate(e11, e12, e, xs11, xs12);
    G4double interpolatedvalue2 = LinLogInterpolate(e21, e22, e, xs21, xs22);
    G4double value = LinLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
  */

  // Lin-Lin
  G4double interpolatedvalue1 = LinLinInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LinLinInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LinLinInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ElasticModel::RandomizeCosTheta(G4double k, const std::size_t& materialID)
{
  G4double integrdiff = 0;  // PROBABILITY between 0 and 1.
  G4double uniformRand = G4UniformRand();
  integrdiff = uniformRand;
  G4double cosTheta = 0.;
  cosTheta = 1 - Theta(G4Electron::ElectronDefinition(), k / eV, integrdiff, materialID);
  // cosTheta = std::cos(theta * CLHEP::pi / 180); ???
  return cosTheta;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNACPA100ElasticModel::ReadDiffCSFile(const std::size_t& materialName,
                                             const G4ParticleDefinition* particleName,
                                             const G4String& file, const G4double&)
{
  const char* path = G4FindDataDir("G4LEDATA");
  if (!path) {
    G4Exception("G4DNACPA100ElasticModel::ReadAllDiffCSFiles", "em0006", FatalException,
                "G4LEDATA environment variable not set.");
    return;
  }

  std::ostringstream fullFileName;
  fullFileName << path << "/" << file << ".dat";

  std::ifstream diffCrossSection(fullFileName.str().c_str());
  // error if file is not there
  std::stringstream endPath;
  if (!diffCrossSection) {
    endPath << "Missing data file: " << file;
    G4Exception("G4DNACPA100ElasticModel::Initialise", "em0003", FatalException,
                endPath.str().c_str());
  }

  tValuesVec[materialName][particleName].push_back(0.);

  G4String line;
  while (std::getline(diffCrossSection, line)) {
    //
    std::istringstream testIss(line);
    G4String test;
    testIss >> test;
    if (test == "#") {
      continue;
    }
    // check if line is empty
    else if (line.empty()) {
      continue;
    }
    std::istringstream iss(line);

    G4double tDummy;
    G4double eDummy;

    iss >> tDummy >> eDummy;

    if (tDummy != tValuesVec[materialName][particleName].back()) {
      // Add the current T value
      tValuesVec[materialName][particleName].push_back(tDummy);
      // Make it correspond to a default zero E value
      eValuesVect[materialName][particleName][tDummy].push_back(0.);
    }
    iss >> diffCrossSectionData[materialName][particleName][tDummy][eDummy];

    if (eDummy != eValuesVect[materialName][particleName][tDummy].back()) {
      eValuesVect[materialName][particleName][tDummy].push_back(eDummy);
    }
  }
}