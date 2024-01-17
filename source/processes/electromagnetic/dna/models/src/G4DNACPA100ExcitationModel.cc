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
// CPA100 excitation model class for electrons
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
// 1/2/2023 : Hoang added modification for DNA cross sections

#include "G4DNACPA100ExcitationModel.hh"

#include "G4DNAChemistryManager.hh"
#include "G4DNAMaterialManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4EnvironmentUtils.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACPA100ExcitationModel::G4DNACPA100ExcitationModel(const G4ParticleDefinition*,
                                                       const G4String& nam)
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNACPA100ExcitationModel::Initialise(const G4ParticleDefinition* p,
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
      G4Exception("G4DNACPA100ExcitationModel::G4DNACPA100ExcitationModel", "CPA001",
                  FatalException, oss.str().c_str());
    }

    const char* path = G4FindDataDir("G4LEDATA");

    if (path == nullptr) {
      G4Exception("G4DNACPA100ExcitationModel::Initialise", "em0006", FatalException,
                  "G4LEDATA environment variable not set.");
      return;
    }

    std::size_t index;
    if (fpG4_WATER != nullptr) {
      index = fpG4_WATER->GetIndex();
      AddCrossSectionData(index, p, "dna/sigma_excitation_e_cpa100", 1.e-20 * m * m);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 255955 * eV);
    }
    if (fpGuanine != nullptr) {
      index = fpGuanine->GetIndex();
      AddCrossSectionData(index, p, "dna/sigma_excitation_e_cpa100_guanine", 1. * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpDeoxyribose != nullptr) {
      index = fpDeoxyribose->GetIndex();
      AddCrossSectionData(index, p, "dna/sigma_excitation_e_cpa100_deoxyribose", 1. * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpCytosine != nullptr) {
      index = fpCytosine->GetIndex();
      AddCrossSectionData(index, p, "dna/sigma_excitation_e_cpa100_cytosine", 1. * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpThymine != nullptr) {
      index = fpThymine->GetIndex();
      AddCrossSectionData(index, p, "dna/sigma_excitation_e_cpa100_thymine", 1. * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpAdenine != nullptr) {
      index = fpAdenine->GetIndex();
      AddCrossSectionData(index, p, "dna/sigma_excitation_e_cpa100_adenine", 1. * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpPhosphate != nullptr) {
      index = fpPhosphate->GetIndex();
      AddCrossSectionData(index, p, "dna/sigma_excitation_e_cpa100_phosphoric_acid", 1. * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }

    LoadCrossSectionData(p);
    G4DNAMaterialManager::Instance()->SetMasterDataModel(DNAModelType::fDNAExcitation, this);
    fpModelData = this;
  }
  else {
    auto dataModel = dynamic_cast<G4DNACPA100ExcitationModel*>(
      G4DNAMaterialManager::Instance()->GetModel(DNAModelType::fDNAExcitation));
    if (dataModel == nullptr) {
      G4cout << "G4DNACPA100ExcitationModel::CrossSectionPerVolume:: not good modelData" << G4endl;
      throw;
    }
    fpModelData = dataModel;
  }
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ExcitationModel::CrossSectionPerVolume(const G4Material* material,
                                                           const G4ParticleDefinition* p,
                                                           G4double ekin, G4double, G4double)
{
  // Get the name of the current particle
  G4String particleName = p->GetParticleName();
  auto MatID = material->GetIndex();
  // initialise variables
  G4double lowLim;
  G4double highLim;
  G4double sigma = 0;

  // Get the low energy limit for the current particle
  lowLim = fpModelData->GetLowELimit(MatID, p);

  // Get the high energy limit for the current particle
  highLim = fpModelData->GetHighELimit(MatID, p);

  // Check that we are in the correct energy range
  if (ekin >= lowLim && ekin < highLim) {
    // Get the map with all the data tables
    auto Data = fpModelData->GetData();

    if ((*Data)[MatID][p] == nullptr) {
      G4Exception("G4DNACPA100ExcitationModel::CrossSectionPerVolume", "em00236", FatalException,
                  "No model is registered");
    }
    // Retrieve the cross section value
    sigma = (*Data)[MatID][p]->FindValue(ekin);

    if (verboseLevel > 2) {
      auto MolDensity =
        (*G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(material))[MatID];
      G4cout << "__________________________________" << G4endl;
      G4cout << "°°° G4DNACPA100ExcitationModel - XS INFO START" << G4endl;
      G4cout << "°°° Kinetic energy(eV)=" << ekin / eV << " particle : " << particleName << G4endl;
      G4cout << "°°° lowLim (eV) = " << lowLim / eV << " highLim (eV) : " << highLim / eV << G4endl;
      G4cout << "°°° Materials = " << (*G4Material::GetMaterialTable())[MatID]->GetName() << G4endl;
      G4cout << "°°° Cross section per " << MatID << " ID molecule (cm^2)=" << sigma / cm / cm
             << G4endl;
      G4cout << "°°° Cross section per Phosphate molecule (cm^-1)="
             << sigma * MolDensity / (1. / cm) << G4endl;
      G4cout << "°°° G4DNACPA100ExcitationModel - XS INFO END" << G4endl;
    }
  }

  // Return the cross section value
  auto MolDensity = (*G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(material))[MatID];
  return sigma * MolDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNACPA100ExcitationModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                                   const G4MaterialCutsCouple* couple,
                                                   const G4DynamicParticle* aDynamicParticle,
                                                   G4double, G4double)
{
  auto materialID = couple->GetMaterial()->GetIndex();
  G4double k = aDynamicParticle->GetKineticEnergy();
  const auto& particle = aDynamicParticle->GetDefinition();
  G4double lowLim = fpModelData->GetLowELimit(materialID, particle);
  G4double highLim = fpModelData->GetHighELimit(materialID, particle);

  // Check if we are in the correct energy range
  if (k >= lowLim && k < highLim) {
    G4int level;
    G4double excitationEnergy;
    G4double newEnergy;
    if (materialID == fpG4_WATER->GetIndex()) {
      level = fpModelData->RandomSelectShell(k, particle, materialID);
      excitationEnergy = eStructure.ExcitationEnergy(level, materialID);
    }
    else {
      do {
        level = eStructure.NumberOfLevels(materialID) * G4UniformRand();
        excitationEnergy = eStructure.ExcitationEnergy(level, materialID);
      } while ((k - eStructure.ExcitationEnergy(level, materialID)) < 0);
    }
    newEnergy = k - excitationEnergy;

    if (k - newEnergy <= 0) {
      G4cout << "k : " << k << "  newEnergy : " << newEnergy << G4endl;
      G4cout << "newEnergy : " << newEnergy << " k : " << k
             << "  excitationEnergy: " << excitationEnergy << G4endl;
      G4cout << "G4DNACPA100ExcitationModel::level : " << eStructure.NumberOfLevels(materialID)
             << " excitationEnergy : " << excitationEnergy << G4endl;
      G4cout << "°°° Materials = " << (*G4Material::GetMaterialTable())[materialID]->GetName()
             << G4endl;
      G4cout << "Attention an error occured !!!" << G4endl;
      abort();
    }

    if (newEnergy >= 0) {
      // We take into account direction change as described page 87 (II.92) in thesis by S. Edel

      G4double cosTheta =
        (excitationEnergy / k) / (1. + (k / (2 * electron_mass_c2)) * (1. - excitationEnergy / k));

      cosTheta = std::sqrt(1. - cosTheta);
      G4double phi = 2. * pi * G4UniformRand();
      const G4ThreeVector& zVers = aDynamicParticle->GetMomentumDirection();
      // Computation of scattering angles (from Subroutine DIRAN in CPA100)

      G4double CT1, ST1, CF1, SF1, CT2, ST2, CF2, SF2;
      G4double sinTheta = std::sqrt(1 - cosTheta * cosTheta);
      CT1 = zVers.z();
      ST1 = std::sqrt(1. - CT1 * CT1);

      ST1 != 0 ? CF1 = zVers.x() / ST1 : CF1 = std::cos(2. * pi * G4UniformRand());
      ST1 != 0 ? SF1 = zVers.y() / ST1 : SF1 = std::sqrt(1. - CF1 * CF1);
      G4double A3, A4, A5, A2, A1;
      A3 = sinTheta * std::cos(phi);
      A4 = A3 * CT1 + ST1 * cosTheta;
      A5 = sinTheta * std::sin(phi);
      A2 = A4 * SF1 + A5 * CF1;
      A1 = A4 * CF1 - A5 * SF1;

      CT2 = CT1 * cosTheta - ST1 * A3;
      ST2 = std::sqrt(1. - CT2 * CT2);

      if (ST2 == 0) {
        ST2 = 1E-6;
      }
      CF2 = A1 / ST2;
      SF2 = A2 / ST2;

      G4ThreeVector zPrimeVers(ST2 * CF2, ST2 * SF2, CT2);
      fParticleChangeForGamma->ProposeMomentumDirection(zPrimeVers.unit());
      if (!statCode) {
        fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
      }
      else {
        fParticleChangeForGamma->SetProposedKineticEnergy(k);
      }

      fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);

      // Chemistry only for water;
      if (materialID == fpG4_WATER->GetIndex()) {
        const G4Track* theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
        G4DNAChemistryManager::Instance()->CreateWaterMolecule(eExcitedMolecule, level,
                                                               theIncomingTrack);
      }
    }
    else {
      G4cerr << "newEnergy : " << newEnergy << " k : " << k
             << "  excitationEnergy: " << excitationEnergy << G4endl;
      G4cerr << "G4DNACPA100ExcitationModel::level : " << eStructure.NumberOfLevels(materialID)
             << " excitationEnergy : " << excitationEnergy << G4endl;
      G4cerr << "°°° Materials = " << (*G4Material::GetMaterialTable())[materialID]->GetName()
             << G4endl;
      G4cerr << "Attention an error occured !!!" << G4endl;
      G4Exception("G4DNACPA100ExcitationModel::SampleSecondaries", "em00236", FatalException,
                  "model is not registered for this energy");
    }
  }
  else {
    G4cerr << "k : " << k << "  lowLim : " << lowLim << "  highLim : " << highLim << G4endl;
    G4Exception("G4DNACPA100ExcitationModel::SampleSecondaries", "em00236", FatalException,
                "model is not registered for this energy");
  }
}
