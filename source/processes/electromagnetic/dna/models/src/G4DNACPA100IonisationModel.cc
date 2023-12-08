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
// CPA100 ionisation model class for electrons
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

#include "G4DNACPA100IonisationModel.hh"

#include "G4DNAChemistryManager.hh"
#include "G4DNAMaterialManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4EnvironmentUtils.hh"
#include "G4LossTableManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACPA100IonisationModel::G4DNACPA100IonisationModel(const G4ParticleDefinition*,
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

void G4DNACPA100IonisationModel::Initialise(const G4ParticleDefinition* p,
                                            const G4DataVector& /*cuts*/)
{
  if (isInitialised) {
    return;
  }
  if (verboseLevel > 3) {
    G4cout << "Calling G4DNACPA100IonisationModel::Initialise()" << G4endl;
  }

  if (!G4DNAMaterialManager::Instance()->IsLocked()) {
    if (p != fpParticle) {
      std::ostringstream oss;
      oss << " Model is not applied for this particle " << p->GetParticleName();
      G4Exception("G4DNACPA100IonisationModel::G4DNACPA100IonisationModel", "CPA001",
                  FatalException, oss.str().c_str());
    }

    const char* path = G4FindDataDir("G4LEDATA");

    if (path == nullptr) {
      G4Exception("G4DNACPA100IonisationModel::Initialise", "em0006", FatalException,
                  "G4LEDATA environment variable not set.");
      return;
    }

    std::size_t index;
    if (fpG4_WATER != nullptr) {
      index = fpG4_WATER->GetIndex();
      G4String eFullFileName = "";
      fasterCode ? eFullFileName = "/dna/sigmadiff_cumulated_ionisation_e_cpa100_rel"
                 : eFullFileName = "/dna/sigmadiff_ionisation_e_cpa100_rel";
      AddCrossSectionData(index, p, "dna/sigma_ionisation_e_cpa100_form_rel", eFullFileName,
                          1.e-20 * m * m);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 255955 * eV);
    }
    if (fpGuanine != nullptr) {
      index = fpGuanine->GetIndex();
      G4String eFullFileName = "";
      if(useDcs) {
        fasterCode ? eFullFileName = "/dna/sigmadiff_cumulated_elastic_e_cpa100_guanine"
                   : eFullFileName = "/dna/sigmadiff_ionisation_e_cpa100_guanine";
      }
      AddCrossSectionData(index, p, "dna/sigma_ionisation_e_cpa100_guanine", eFullFileName,
                          1. * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpDeoxyribose != nullptr) {
      index = fpDeoxyribose->GetIndex();
      G4String eFullFileName = "";
      if(useDcs) {
        eFullFileName = "/dna/sigmadiff_cumulated_ionisation_e_cpa100_deoxyribose";
      }
      AddCrossSectionData(index, p, "dna/sigma_ionisation_e_cpa100_deoxyribose", eFullFileName,
                          1. * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpCytosine != nullptr) {
      index = fpCytosine->GetIndex();
      G4String eFullFileName = "";
      if(useDcs) {
        fasterCode ? eFullFileName = "/dna/sigmadiff_cumulated_ionisation_e_cpa100_cytosine"
                   : eFullFileName = "/dna/sigmadiff_ionisation_e_cpa100_cytosine";
      }
      AddCrossSectionData(index, p, "dna/sigma_ionisation_e_cpa100_cytosine", eFullFileName,
                          1. * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpThymine != nullptr) {
      index = fpThymine->GetIndex();
      G4String eFullFileName = "";
      if(useDcs) {
        fasterCode ? eFullFileName = "/dna/sigmadiff_cumulated_ionisation_e_cpa100_thymine"
                   : eFullFileName = "/dna/sigmadiff_ionisation_e_cpa100_thymine";
      }
      AddCrossSectionData(index, p, "dna/sigma_ionisation_e_cpa100_thymine", eFullFileName,
                          1. * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpAdenine != nullptr) {
      index = fpAdenine->GetIndex();
      G4String eFullFileName = "";
      if(useDcs) {
        fasterCode ? eFullFileName = "/dna/sigmadiff_cumulated_ionisation_e_cpa100_adenine"
                   : eFullFileName = "/dna/sigmadiff_ionisation_e_cpa100_adenine";
      }
      AddCrossSectionData(index, p, "dna/sigma_ionisation_e_cpa100_adenine", eFullFileName,
                          1. * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    if (fpPhosphate != nullptr) {
      index = fpPhosphate->GetIndex();
      G4String eFullFileName = "";
      if(useDcs) {
        eFullFileName = "dna/sigmadiff_cumulated_ionisation_e_cpa100_phosphoric_acid";
      }
      AddCrossSectionData(index, p, "dna/sigma_ionisation_e_cpa100_phosphoric_acid",eFullFileName,
                          1. * cm * cm);
      SetLowELimit(index, p, 11 * eV);
      SetHighELimit(index, p, 1 * MeV);
    }
    LoadCrossSectionData(p);
    G4DNAMaterialManager::Instance()->SetMasterDataModel(DNAModelType::fDNAIonisation, this);
    fpModelData = this;
  }
  else {
    auto dataModel = dynamic_cast<G4DNACPA100IonisationModel*>(
      G4DNAMaterialManager::Instance()->GetModel(DNAModelType::fDNAIonisation));
    if (dataModel == nullptr) {
      G4cout << "G4DNACPA100IonisationModel::CrossSectionPerVolume:: not good modelData" << G4endl;
      throw;
    }
    fpModelData = dataModel;
  }

  fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();

  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

G4double G4DNACPA100IonisationModel::CrossSectionPerVolume(const G4Material* material,
                                                           const G4ParticleDefinition* p,
                                                           G4double ekin, G4double, G4double)
{
  // initialise the cross section value (output value)
  G4double sigma(0);

  // Get the current particle name
  const G4String& particleName = p->GetParticleName();

  if (p != fpParticle) {
    G4Exception("G4DNACPA100IonisationModel::CrossSectionPerVolume", "em00223", FatalException,
                "No model is registered for this particle");
  }

  auto matID = material->GetIndex();

  // Set the low and high energy limits
  G4double lowLim = fpModelData->GetLowELimit(matID, p);
  G4double highLim = fpModelData->GetHighELimit(matID, p);

  // Check that we are in the correct energy range
  if (ekin >= lowLim && ekin < highLim) {
    // Get the map with all the model data tables
    auto tableData = fpModelData->GetData();

    if ((*tableData)[matID][p] == nullptr) {
      G4Exception("G4DNACPA100IonisationModel::CrossSectionPerVolume", "em00236", FatalException,
                  "No model is registered");
    }
    else {
      sigma = (*tableData)[matID][p]->FindValue(ekin);
    }

    if (verboseLevel > 2) {
      auto MolDensity =
        (*G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(material))[matID];
      G4cout << "__________________________________" << G4endl;
      G4cout << "°°° G4DNACPA100IonisationModel - XS INFO START" << G4endl;
      G4cout << "°°° Kinetic energy(eV)=" << ekin / eV << " particle : " << particleName << G4endl;
      G4cout << "°°° lowLim (eV) = " << lowLim / eV << " highLim (eV) : " << highLim / eV << G4endl;
      G4cout << "°°° Materials = " << (*G4Material::GetMaterialTable())[matID]->GetName() << G4endl;
      G4cout << "°°° Cross section per " << matID << " index molecule (cm^2)=" << sigma / cm / cm
             << G4endl;
      G4cout << "°°° Cross section per Phosphate molecule (cm^-1)="
             << sigma * MolDensity / (1. / cm) << G4endl;
      G4cout << "°°° G4DNACPA100IonisationModel - XS INFO END" << G4endl;
    }
  }

  auto MolDensity = (*G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(material))[matID];
  return sigma * MolDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNACPA100IonisationModel::SampleSecondaries(
  std::vector<G4DynamicParticle*>* fvect,
  const G4MaterialCutsCouple* couple,  // must be set!
  const G4DynamicParticle* particle, G4double, G4double)
{
  if (verboseLevel > 3) {
    G4cout << "Calling SampleSecondaries() of G4DNACPA100IonisationModel" << G4endl;
  }
  auto k = particle->GetKineticEnergy();

  const G4Material* material = couple->GetMaterial();

  auto MatID = material->GetIndex();

  auto p = particle->GetDefinition();

  auto lowLim = fpModelData->GetLowELimit(MatID, p);
  auto highLim = fpModelData->GetHighELimit(MatID, p);

  // Check if we are in the correct energy range
  if (k >= lowLim && k < highLim) {
    const auto& primaryDirection = particle->GetMomentumDirection();
    auto particleMass = particle->GetDefinition()->GetPDGMass();
    auto totalEnergy = k + particleMass;
    auto pSquare = k * (totalEnergy + particleMass);
    auto totalMomentum = std::sqrt(pSquare);
    G4int shell = -1;
    G4double bindingEnergy, secondaryKinetic;
    shell = fpModelData->RandomSelectShell(k, p, MatID);
    bindingEnergy = iStructure.IonisationEnergy(shell, MatID);

    if (k < bindingEnergy) {
      return;
    }

    auto info = std::make_tuple(MatID, k, shell);

    secondaryKinetic = -1000 * eV;
    if (fpG4_WATER->GetIndex() != MatID) {//for DNA material useDcs = false
      secondaryKinetic = fpModelData->RandomizeEjectedElectronEnergyFromanalytical(info);
    }else if(fasterCode){
        secondaryKinetic = fpModelData->RandomizeEjectedElectronEnergyFromCumulatedDcs(info);
      }else{
        secondaryKinetic = fpModelData->RandomizeEjectedElectronEnergy(info);
      }

    G4double cosTheta = 0.;
    G4double phi = 0.;
    RandomizeEjectedElectronDirection(particle->GetDefinition(), k, secondaryKinetic, cosTheta,
                                      phi);

    G4double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
    G4double dirX = sinTheta * std::cos(phi);
    G4double dirY = sinTheta * std::sin(phi);
    G4double dirZ = cosTheta;
    G4ThreeVector deltaDirection(dirX, dirY, dirZ);
    deltaDirection.rotateUz(primaryDirection);

    // SI - For atom. deexc. tagging - 23/05/2017
    if (secondaryKinetic > 0) {
      auto dp = new G4DynamicParticle(G4Electron::Electron(), deltaDirection, secondaryKinetic);
      fvect->push_back(dp);
    }

    if (particle->GetDefinition() != fpParticle) {
      fParticleChangeForGamma->ProposeMomentumDirection(primaryDirection);
    }
    else {
      G4double deltaTotalMomentum =
        std::sqrt(secondaryKinetic * (secondaryKinetic + 2. * electron_mass_c2));
      G4double finalPx =
        totalMomentum * primaryDirection.x() - deltaTotalMomentum * deltaDirection.x();
      G4double finalPy =
        totalMomentum * primaryDirection.y() - deltaTotalMomentum * deltaDirection.y();
      G4double finalPz =
        totalMomentum * primaryDirection.z() - deltaTotalMomentum * deltaDirection.z();
      G4double finalMomentum = std::sqrt(finalPx * finalPx + finalPy * finalPy + finalPz * finalPz);
      finalPx /= finalMomentum;
      finalPy /= finalMomentum;
      finalPz /= finalMomentum;

      G4ThreeVector direction;
      direction.set(finalPx, finalPy, finalPz);

      fParticleChangeForGamma->ProposeMomentumDirection(direction.unit());
    }

    // SI - For atom. deexc. tagging - 23/05/2017

    // AM: sample deexcitation
    // here we assume that H_{2}O electronic levels are the same of Oxigen.
    // this can be considered true with a rough 10% error in energy on K-shell,

    G4double scatteredEnergy = k - bindingEnergy - secondaryKinetic;

    // SI: only atomic deexcitation from K shell is considered
    // Hoang: only for water
    if (fpG4_WATER != nullptr && material == G4Material::GetMaterial("G4_WATER")) {
      std::size_t secNumberInit = 0;  // need to know at a certain point the energy of secondaries
      std::size_t secNumberFinal = 0;  // So I'll make the diference and then sum the energies
      if ((fAtomDeexcitation != nullptr) && shell == 4) {
        G4int Z = 8;
        auto Kshell = fAtomDeexcitation->GetAtomicShell(Z, G4AtomicShellEnumerator(0));
        secNumberInit = fvect->size();
        fAtomDeexcitation->GenerateParticles(fvect, Kshell, Z, 0, 0);
        secNumberFinal = fvect->size();
        if (secNumberFinal > secNumberInit) {
          for (std::size_t i = secNumberInit; i < secNumberFinal; ++i) {
            // Check if there is enough residual energy
            if (bindingEnergy >= ((*fvect)[i])->GetKineticEnergy()) {
              // Ok, this is a valid secondary: keep it
              bindingEnergy -= ((*fvect)[i])->GetKineticEnergy();
            }
            else {
              // Invalid secondary: not enough energy to create it!
              // Keep its energy in the local deposit
              delete (*fvect)[i];
              (*fvect)[i] = nullptr;
            }
          }
        }
      }
    }

    // This should never happen
    if (bindingEnergy < 0.0) {
      G4Exception("G4DNACPA100IonisatioModel1::SampleSecondaries()", "em2050", FatalException,
                  "Negative local energy deposit");
    }
    if (!statCode) {
      fParticleChangeForGamma->SetProposedKineticEnergy(scatteredEnergy);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(bindingEnergy);
    }
    else {
      fParticleChangeForGamma->SetProposedKineticEnergy(k);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(k - scatteredEnergy);
    }

    // only water for chemistry
    if (fpG4_WATER != nullptr && material == G4Material::GetMaterial("G4_WATER")) {
      const G4Track* theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
      G4DNAChemistryManager::Instance()->CreateWaterMolecule(eIonizedMolecule, shell,
                                                             theIncomingTrack);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNACPA100IonisationModel::RandomizeEjectedElectronEnergy(PartKineticInMat info)
{
  auto MatID = std::get<0>(info);
  auto k = std::get<1>(info);
  auto shell = std::get<2>(info);
  G4double maximumEnergyTransfer = 0.;
  auto IonLevel = iStructure.IonisationEnergy(shell, MatID);
  (k + IonLevel) / 2. > k ? maximumEnergyTransfer = k : maximumEnergyTransfer = (k + IonLevel) / 2.;

  G4double crossSectionMaximum = 0.;

  G4double minEnergy = IonLevel;
  G4double maxEnergy = maximumEnergyTransfer;

  // nEnergySteps can be optimized - 100 by default
  G4int nEnergySteps = 50;

  G4double value(minEnergy);
  G4double stpEnergy(std::pow(maxEnergy / value, 1. / static_cast<G4double>(nEnergySteps - 1)));
  G4int step(nEnergySteps);
  G4double differentialCrossSection = 0.;
  while (step > 0) {
    step--;
    differentialCrossSection = DifferentialCrossSection(info, value / eV);

    if (differentialCrossSection > 0) {
      crossSectionMaximum = differentialCrossSection;
      break;
    }
    value *= stpEnergy;
  }

  G4double secondaryElectronKineticEnergy = 0.;
  do {
    secondaryElectronKineticEnergy = G4UniformRand() * (maximumEnergyTransfer - IonLevel);
  } while (G4UniformRand() * crossSectionMaximum
           > DifferentialCrossSection(info, (secondaryElectronKineticEnergy + IonLevel) / eV));

  return secondaryElectronKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNACPA100IonisationModel::RandomizeEjectedElectronDirection(G4ParticleDefinition*,
                                                                   G4double k, G4double secKinetic,
                                                                   G4double& cosTheta,
                                                                   G4double& phi)
{
  phi = twopi * G4UniformRand();
  G4double sin2O = (1. - secKinetic / k) / (1. + secKinetic / (2. * electron_mass_c2));
  cosTheta = std::sqrt(1. - sin2O);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNACPA100IonisationModel::DifferentialCrossSection(PartKineticInMat info,
                                                              const G4double& energyTransfer)
{
  auto MatID = std::get<0>(info);
  auto k = std::get<1>(info) / eV;  // in eV unit
  auto shell = std::get<2>(info);
  G4double sigma = 0.;
  G4double shellEnergy = iStructure.IonisationEnergy(shell, MatID);
  G4double kSE(energyTransfer - shellEnergy);

  if (energyTransfer >= shellEnergy) {
    G4double valueT1 = 0;
    G4double valueT2 = 0;
    G4double valueE21 = 0;
    G4double valueE22 = 0;
    G4double valueE12 = 0;
    G4double valueE11 = 0;

    G4double xs11 = 0;
    G4double xs12 = 0;
    G4double xs21 = 0;
    G4double xs22 = 0;

    auto t2 = std::upper_bound(fTMapWithVec[MatID][fpParticle].begin(),
                               fTMapWithVec[MatID][fpParticle].end(), k);
    auto t1 = t2 - 1;

    if (kSE <= fEMapWithVector[MatID][fpParticle][(*t1)].back()
        && kSE <= fEMapWithVector[MatID][fpParticle][(*t2)].back())
    {
      auto e12 = std::upper_bound(fEMapWithVector[MatID][fpParticle][(*t1)].begin(),
                                  fEMapWithVector[MatID][fpParticle][(*t1)].end(), kSE);
      auto e11 = e12 - 1;

      auto e22 = std::upper_bound(fEMapWithVector[MatID][fpParticle][(*t2)].begin(),
                                  fEMapWithVector[MatID][fpParticle][(*t2)].end(), kSE);
      auto e21 = e22 - 1;

      valueT1 = *t1;
      valueT2 = *t2;
      valueE21 = *e21;
      valueE22 = *e22;
      valueE12 = *e12;
      valueE11 = *e11;

      xs11 = diffCrossSectionData[MatID][fpParticle][shell][valueT1][valueE11];
      xs12 = diffCrossSectionData[MatID][fpParticle][shell][valueT1][valueE12];
      xs21 = diffCrossSectionData[MatID][fpParticle][shell][valueT2][valueE21];
      xs22 = diffCrossSectionData[MatID][fpParticle][shell][valueT2][valueE22];
    }

    G4double xsProduct = xs11 * xs12 * xs21 * xs22;

    if (xsProduct != 0.) {
      sigma = QuadInterpolator(valueE11, valueE12, valueE21, valueE22, xs11, xs12, xs21, xs22,
                               valueT1, valueT2, k, kSE);
    }
  }

  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNACPA100IonisationModel::Interpolate(G4double e1, G4double e2, G4double e, G4double xs1,
                                                 G4double xs2)
{
  G4double value = 0.;

  // Log-log interpolation by default

  if (e1 != 0 && e2 != 0 && (std::log10(e2) - std::log10(e1)) != 0 && !fasterCode) {
    G4double a = (std::log10(xs2) - std::log10(xs1)) / (std::log10(e2) - std::log10(e1));
    G4double b = std::log10(xs2) - a * std::log10(e2);
    G4double sigma = a * std::log10(e) + b;
    value = (std::pow(10., sigma));
  }

  // Switch to lin-lin interpolation
  /*
  if ((e2-e1)!=0)
  {
    G4double d1 = xs1;
    G4double d2 = xs2;
    value = (d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
  }
  */

  // Switch to log-lin interpolation for faster code

  if ((e2 - e1) != 0 && xs1 != 0 && xs2 != 0 && fasterCode) {
    G4double d1 = std::log10(xs1);
    G4double d2 = std::log10(xs2);
    value = std::pow(10., (d1 + (d2 - d1) * (e - e1) / (e2 - e1)));
  }

  // Switch to lin-lin interpolation for faster code
  // in case one of xs1 or xs2 (=cum proba) value is zero

  if ((e2 - e1) != 0 && (xs1 == 0 || xs2 == 0) && fasterCode) {
    G4double d1 = xs1;
    G4double d2 = xs2;
    value = (d1 + (d2 - d1) * (e - e1) / (e2 - e1));
  }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNACPA100IonisationModel::QuadInterpolator(G4double e11, G4double e12, G4double e21,
                                                      G4double e22, G4double xs11, G4double xs12,
                                                      G4double xs21, G4double xs22, G4double t1,
                                                      G4double t2, G4double t, G4double e)
{
  G4double interpolatedvalue1 = Interpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = Interpolate(e21, e22, e, xs21, xs22);
  G4double value = Interpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);

  return value;
}

G4double
G4DNACPA100IonisationModel::RandomizeEjectedElectronEnergyFromCumulatedDcs(PartKineticInMat info)
{
  auto MatID = std::get<0>(info);
  auto shell = std::get<2>(info);
  G4double secondaryElectronKineticEnergy =
    RandomTransferedEnergy(info) * eV - iStructure.IonisationEnergy(shell, MatID);
  if (secondaryElectronKineticEnergy < 0.) {
    return 0.;
  }
  return secondaryElectronKineticEnergy;
}

G4double G4DNACPA100IonisationModel::RandomTransferedEnergy(PartKineticInMat info)
{
  auto materialID = std::get<0>(info);
  auto k = std::get<1>(info) / eV;  // data table in eV
  auto shell = std::get<2>(info);
  G4double ejectedElectronEnergy = 0.;
  G4double valueK1 = 0;
  G4double valueK2 = 0;
  G4double valueCumulCS21 = 0;
  G4double valueCumulCS22 = 0;
  G4double valueCumulCS12 = 0;
  G4double valueCumulCS11 = 0;
  G4double secElecE11 = 0;
  G4double secElecE12 = 0;
  G4double secElecE21 = 0;
  G4double secElecE22 = 0;

  if (k == fTMapWithVec[materialID][fpParticle].back()) {
    k = k * (1. - 1e-12);
  }

  G4double random = G4UniformRand();
  auto k2 = std::upper_bound(fTMapWithVec[materialID][fpParticle].begin(),
                             fTMapWithVec[materialID][fpParticle].end(), k);
  auto k1 = k2 - 1;

  if (random <= fProbaShellMap[materialID][fpParticle][shell][(*k1)].back()
      && random <= fProbaShellMap[materialID][fpParticle][shell][(*k2)].back())
  {
    auto cumulCS12 =
      std::upper_bound(fProbaShellMap[materialID][fpParticle][shell][(*k1)].begin(),
                       fProbaShellMap[materialID][fpParticle][shell][(*k1)].end(), random);
    auto cumulCS11 = cumulCS12 - 1;
    // Second one.
    auto cumulCS22 =
      std::upper_bound(fProbaShellMap[materialID][fpParticle][shell][(*k2)].begin(),
                       fProbaShellMap[materialID][fpParticle][shell][(*k2)].end(), random);
    auto cumulCS21 = cumulCS22 - 1;

    valueK1 = *k1;
    valueK2 = *k2;
    valueCumulCS11 = *cumulCS11;
    valueCumulCS12 = *cumulCS12;
    valueCumulCS21 = *cumulCS21;
    valueCumulCS22 = *cumulCS22;

    secElecE11 = fEnergySecondaryData[materialID][fpParticle][shell][valueK1][valueCumulCS11];
    secElecE12 = fEnergySecondaryData[materialID][fpParticle][shell][valueK1][valueCumulCS12];
    secElecE21 = fEnergySecondaryData[materialID][fpParticle][shell][valueK2][valueCumulCS21];
    secElecE22 = fEnergySecondaryData[materialID][fpParticle][shell][valueK2][valueCumulCS22];

    if (valueCumulCS11 == 0. && valueCumulCS12 == 1.) {
      auto interpolatedvalue2 =
        Interpolate(valueCumulCS21, valueCumulCS22, random, secElecE21, secElecE22);
      G4double valueNrjTransf = Interpolate(valueK1, valueK2, k, 0., interpolatedvalue2);
      return valueNrjTransf;
    }
  }

  if (random > fProbaShellMap[materialID][fpParticle][shell][(*k1)].back()) {
    auto cumulCS22 =
      std::upper_bound(fProbaShellMap[materialID][fpParticle][shell][(*k2)].begin(),
                       fProbaShellMap[materialID][fpParticle][shell][(*k2)].end(), random);
    auto cumulCS21 = cumulCS22 - 1;
    valueK1 = *k1;
    valueK2 = *k2;
    valueCumulCS21 = *cumulCS21;
    valueCumulCS22 = *cumulCS22;

    secElecE21 = fEnergySecondaryData[materialID][fpParticle][shell][valueK2][valueCumulCS21];
    secElecE22 = fEnergySecondaryData[materialID][fpParticle][shell][valueK2][valueCumulCS22];

    G4double interpolatedvalue2 =
      Interpolate(valueCumulCS21, valueCumulCS22, random, secElecE21, secElecE22);

    G4double value = Interpolate(valueK1, valueK2, k, 0., interpolatedvalue2);
    return value;
  }
  G4double nrjTransfProduct = secElecE11 * secElecE12 * secElecE21 * secElecE22;

  if (nrjTransfProduct != 0.) {
    ejectedElectronEnergy =
      QuadInterpolator(valueCumulCS11, valueCumulCS12, valueCumulCS21, valueCumulCS22, secElecE11,
                       secElecE12, secElecE21, secElecE22, valueK1, valueK2, k, random);
  }
  return ejectedElectronEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4DNACPA100IonisationModel::RandomizeEjectedElectronEnergyFromanalytical(PartKineticInMat info)
{
  auto MatID = std::get<0>(info);
  auto tt = std::get<1>(info);
  auto shell = std::get<2>(info);
  // ***** METHOD by M. C. Bordage ***** (optimized)
  //  Composition sampling method based on eq 7 in (Guerra et al. 2015) the RBEBV

  //// Defining constants
  G4double alfa = 1. / 137;  // fine structure constant
  G4double e_charge = 1.6e-19;  // electron charge
  G4double e_mass = 9.1e-31;  // electron mass in kg
  G4double c = 3e8;  // speed of light in vacuum constant c (m/s)
  G4double mc2 = e_mass * c * c / e_charge;  //

  G4double BB = iStructure.IonisationEnergy(shell, MatID);  // binding energy of the shell (eV)

  if (tt <= BB) return 0.;

  G4double b_prime = BB / mc2;  // binding energy divided by mc2
  G4double beta_b2 = 1. - 1. / ((1 + b_prime) * (1 + b_prime));  // binding energy Beta

  //// Indicent energy
  //// tt is the incident electron energy

  G4double t_prime = tt / mc2;  // incident energy divided by mc2
  G4double t = tt / BB;  // reduced incident energy by binding energy

  G4double D = (1 + 2 * t_prime) / ((1 + t_prime / 2) * (1 + t_prime / 2));
  G4double F = b_prime * b_prime / ((1 + t_prime / 2) * (1 + t_prime / 2));

  G4double beta_t2 = 1 - 1 / ((1 + t_prime) * (1 + t_prime));  // incident energy Beta

  G4double PHI_R = std::cos(std::sqrt(alfa * alfa / (beta_t2 + beta_b2))
                            * std::log(beta_t2 / beta_b2));  // relativistic Vriens function phi
  G4double G_R = std::log(beta_t2 / (1 - beta_t2)) - beta_t2 - std::log(2 * b_prime);

  G4double tplus1 = t + 1;
  G4double tminus1 = t - 1;
  G4double tplus12 = tplus1 * tplus1;
  G4double ZH1max = 1 + F - (PHI_R * D * (2 * t + 1) / (2 * t * tplus1));
  G4double ZH2max = 1 - PHI_R * D / 4;

  G4double A1_p = ZH1max * tminus1 / tplus1;  // A1'
  G4double A2_p = ZH2max * tminus1 / (t * tplus1);  // A2'
  G4double A3_p = ((tplus12 - 4) / tplus12) * G_R;  // A3'

  G4double AAA = A1_p + A2_p + A3_p;

  G4double AA1_R = A1_p / AAA;
  G4double AA2_R = (A1_p + A2_p) / AAA;

  G4int FF = 0;
  G4double fx = 0;
  G4double gx = 0;
  G4double gg = 0;
  G4double wx = 0;

  G4double r1 = 0;
  G4double r2 = 0;
  G4double r3 = 0;

  //

  do {
    r1 = G4UniformRand();
    r2 = G4UniformRand();
    r3 = G4UniformRand();

    if (r1 > AA2_R)
      FF = 3;
    else if ((r1 > AA1_R) && (r1 < AA2_R))
      FF = 2;
    else
      FF = 1;

    switch (FF) {
      case 1: {
        fx = r2 * tminus1 / tplus1;
        wx = 1 / (1 - fx) - 1;
        gg = PHI_R * D * (wx + 1) / tplus1;
        gx = 1 - gg;
        gx = gx - gg * (wx + 1) / (2 * (t - wx));
        gx = gx + F * (wx + 1) * (wx + 1);
        gx = gx / ZH1max;
        break;
      }

      case 2: {
        fx = tplus1 + r2 * tminus1;
        wx = t * tminus1 * r2 / fx;
        gx = 1 - (PHI_R * D * (t - wx) / (2 * tplus1));
        gx = gx / ZH2max;
        break;
      }

      case 3: {
        fx = 1 - r2 * (tplus12 - 4) / tplus12;
        wx = std::sqrt(1 / fx) - 1;
        gg = (wx + 1) / (t - wx);
        gx = (1 + gg * gg * gg) / 2;
        break;
      }
    }  // switch

  } while (r3 > gx);

  return wx * BB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNACPA100IonisationModel::ReadDiffCSFile(const std::size_t& materialID,
                                                const G4ParticleDefinition* p, const G4String& file,
                                                const G4double& scaleFactor)
{
  const char* path = G4FindDataDir("G4LEDATA");
  if (path == nullptr) {
    G4Exception("G4DNACPA100IonisationModel::ReadAllDiffCSFiles", "em0006", FatalException,
                "G4LEDATA environment variable was not set.");
    return;
  }

  std::ostringstream fullFileName;
  fullFileName << path << "/" << file << ".dat";

  std::ifstream diffCrossSection(fullFileName.str().c_str());
  std::stringstream endPath;
  if (!diffCrossSection) {
    endPath << "Missing data file: " << file;
    G4Exception("G4DNACPA100IonisationModel::Initialise", "em0003", FatalException,
                endPath.str().c_str());
  }

  // load data from the file
  fTMapWithVec[materialID][p].push_back(0.);

  G4String line;

  while (!diffCrossSection.eof()) {
    G4double T, E;
    diffCrossSection >> T >> E;

    if (T != fTMapWithVec[materialID][p].back()) {
      fTMapWithVec[materialID][p].push_back(T);
    }

    // T is incident energy, E is the energy transferred
    if (T != fTMapWithVec[materialID][p].back()) {
      fTMapWithVec[materialID][p].push_back(T);
    }

    auto eshell = (G4int)iStructure.NumberOfLevels(materialID);
    for (G4int shell = 0; shell < eshell; ++shell) {
      diffCrossSection >> diffCrossSectionData[materialID][p][shell][T][E];
      if (fasterCode) {
        fEnergySecondaryData[materialID][p][shell][T]
                            [diffCrossSectionData[materialID][p][shell][T][E]] = E;

        fProbaShellMap[materialID][p][shell][T].push_back(
          diffCrossSectionData[materialID][p][shell][T][E]);
      }
      else {
        diffCrossSectionData[materialID][p][shell][T][E] *= scaleFactor;
        fEMapWithVector[materialID][p][T].push_back(E);
      }
    }
  }
}
