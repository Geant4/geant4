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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Models come from
// M. Bug et al, Rad. Phys and Chem. 130, 459-479 (2017)
//

#include "G4DNAPTBExcitationModel.hh"

#include "G4DNAChemistryManager.hh"
#include "G4DNAMaterialManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4SystemOfUnits.hh"

G4DNAPTBExcitationModel::G4DNAPTBExcitationModel(const G4String& applyToMaterial,
                                                 const G4ParticleDefinition*, const G4String& nam)
  : G4VDNAModel(nam, applyToMaterial)
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
  // initialisation of mean energy loss for each material

  if (fpTHF != nullptr) {
    fTableMeanEnergyPTB[fpTHF->GetIndex()] = 8.01 * eV;
  }

  if (fpPY != nullptr) {
    fTableMeanEnergyPTB[fpPY->GetIndex()] = 7.61 * eV;
  }

  if (fpPU != nullptr) {
    fTableMeanEnergyPTB[fpPU->GetIndex()] = 7.61 * eV;
  }
  if (fpTMP != nullptr) {
    fTableMeanEnergyPTB[fpTMP->GetIndex()] = 8.01 * eV;
  }

  if (verboseLevel > 0) {
    G4cout << "PTB excitation model is constructed " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBExcitationModel::Initialise(const G4ParticleDefinition* particle,
                                         const G4DataVector& /*cuts*/)
{
  if (isInitialised) {
    return;
  }
  if (verboseLevel > 3)
  {
    G4cout << "Calling G4DNAPTBExcitationModel::Initialise()" << G4endl;
  }

  if (particle != G4Electron::ElectronDefinition()) {
    std::ostringstream oss;
    oss << " Model is not applied for this particle " << particle->GetParticleName();
    G4Exception("G4DNAPTBExcitationModel::Initialise", "PTB001", FatalException, oss.str().c_str());
  }

  G4double scaleFactor = 1e-16 * cm * cm;
  G4double scaleFactorBorn = (1.e-22 / 3.343) * m * m;

  //*******************************************************
  // Cross section data
  //*******************************************************
  std::size_t index;
  if (fpTHF != nullptr) {
    index = fpTHF->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_excitation_e-_PTB_THF", scaleFactor);
    SetLowELimit(index, particle, 9. * eV);
    SetHighELimit(index, particle, 1. * keV);
  }
  if (fpPY != nullptr) {
    index = fpPY->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_excitation_e-_PTB_PY", scaleFactor);
    SetLowELimit(index, particle, 9. * eV);
    SetHighELimit(index, particle, 1. * keV);
  }

  if (fpPU != nullptr) {
    index = fpPU->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_excitation_e-_PTB_PU", scaleFactor);
    SetLowELimit(index, particle, 9. * eV);
    SetHighELimit(index, particle, 1. * keV);
  }

  if (fpTMP != nullptr) {
    index = fpTMP->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_excitation_e-_PTB_TMP", scaleFactor);
    SetLowELimit(index, particle, 9. * eV);
    SetHighELimit(index, particle, 1. * keV);
  }
  if (fpG4_WATER != nullptr) {
    index = fpG4_WATER->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_excitation_e_born", scaleFactorBorn);
    SetLowELimit(index, particle, 9. * eV);
    SetHighELimit(index, particle, 1. * keV);
  }
  // DNA materials
  //
  if (fpBackbone_THF != nullptr) {
    index = fpBackbone_THF->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_excitation_e-_PTB_THF", scaleFactor * 33. / 30);
    SetLowELimit(index, particle, 9. * eV);
    SetHighELimit(index, particle, 1. * keV);
  }
  if (fpCytosine_PY != nullptr) {
    index = fpCytosine_PY->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_excitation_e-_PTB_PY", scaleFactor * 42. / 30);
    SetLowELimit(index, particle, 9. * eV);
    SetHighELimit(index, particle, 1. * keV);
  }
  if (fpThymine_PY != nullptr) {
    index = fpThymine_PY->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_excitation_e-_PTB_PY", scaleFactor * 48. / 30);
    SetLowELimit(index, particle, 9. * eV);
    SetHighELimit(index, particle, 1. * keV);
  }
  if (fpAdenine_PU != nullptr) {
    index = fpAdenine_PU->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_excitation_e-_PTB_PU", scaleFactor * 50. / 44);
    SetLowELimit(index, particle, 9. * eV);
    SetHighELimit(index, particle, 1. * keV);
  }
  if (fpGuanine_PU != nullptr) {
    index = fpGuanine_PU->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_excitation_e-_PTB_PU", scaleFactor * 56. / 44);
    SetLowELimit(index, particle, 9. * eV);
    SetHighELimit(index, particle, 1. * keV);
  }
  if (fpBackbone_TMP != nullptr) {
    index = fpBackbone_TMP->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_excitation_e-_PTB_TMP", scaleFactor * 33. / 50);
    SetLowELimit(index, particle, 9. * eV);
    SetHighELimit(index, particle, 1. * keV);
  }
  // MPietrzak, adding paths for N2
  if (fpN2 != nullptr) {
    index = fpN2->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_excitation_e-_PTB_N2", scaleFactor);
    SetLowELimit(index, particle, 13. * eV);
    SetHighELimit(index, particle, 1.02 * MeV);
  }
  if (!G4DNAMaterialManager::Instance()->IsLocked()) {
    // Load data
    LoadCrossSectionData(particle);
    G4DNAMaterialManager::Instance()->SetMasterDataModel(DNAModelType::fDNAExcitation, this);
    fpModelData = this;
  }
  else {
    auto dataModel = dynamic_cast<G4DNAPTBExcitationModel*>(
      G4DNAMaterialManager::Instance()->GetModel(DNAModelType::fDNAExcitation));
    if (dataModel == nullptr) {
      G4cout << "G4DNAPTBExcitationModel::Initialise:: not good modelData" << G4endl;
      G4Exception("G4DNAPTBExcitationModel::Initialise", "PTB0006", FatalException,
                  "not good modelData");
    }
    else {
      fpModelData = dataModel;
    }
  }

  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBExcitationModel::CrossSectionPerVolume(const G4Material* material,
                                                        const G4ParticleDefinition* p,
                                                        G4double ekin, G4double /*emin*/,
                                                        G4double /*emax*/)
{
  // Get the name of the current particle
  G4String particleName = p->GetParticleName();
  const std::size_t& MatID = material->GetIndex();
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
      G4Exception("G4DNAPTBExcitationModel::CrossSectionPerVolume", "em00236", FatalException,
                  "No model is registered");
    }
    // Retrieve the cross section value
    sigma = (*Data)[MatID][p]->FindValue(ekin);

    if (verboseLevel > 2) {
      G4cout << "__________________________________" << G4endl;
      G4cout << "°°° G4DNAPTBExcitationModel - XS INFO START" << G4endl;
      G4cout << "°°° Kinetic energy(eV)=" << ekin / eV << " particle : " << particleName << G4endl;
      G4cout << "°°° Cross section per " << MatID << " ID molecule (cm^2)=" << sigma / cm / cm
             << G4endl;
      G4cout << "°°° G4DNAPTBExcitationModel - XS INFO END" << G4endl;
    }
  }

  // Return the cross section value
  auto MolDensity =
    (*G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(material))[MatID];
  return sigma * MolDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBExcitationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                const G4MaterialCutsCouple* couple,
                                                const G4DynamicParticle* aDynamicParticle,
                                                G4double /*tmin*/, G4double /*tmax*/)
{
  const std::size_t& materialID = (std::size_t)couple->GetIndex();

  // Get the incident particle kinetic energy
  G4double k = aDynamicParticle->GetKineticEnergy();
  // Get the particle name
  const auto& particle = aDynamicParticle->GetDefinition();
  // Get the energy limits
  G4double lowLim = fpModelData->GetLowELimit(materialID, particle);
  G4double highLim = fpModelData->GetHighELimit(materialID, particle);

  // Check if we are in the correct energy range
  if (k >= lowLim && k < highLim) {
    if (fpN2 != nullptr && materialID == fpN2->GetIndex()) {
      // Retrieve the excitation energy for the current material
      G4int level = fpModelData->RandomSelectShell(k, particle, materialID);
      G4double excitationEnergy = ptbExcitationStructure.ExcitationEnergy(level, fpN2->GetIndex());

      // Calculate the new energy of the particle
      G4double newEnergy = k - excitationEnergy;

      // Check that the new energy is above zero before applying it the particle.
      // Otherwise, do nothing.
      if (newEnergy > 0) {
        fParticleChangeForGamma->ProposeMomentumDirection(aDynamicParticle->GetMomentumDirection());
        fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
        fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
        G4double ioniThres = ptbIonisationStructure.IonisationEnergy(0, fpN2->GetIndex());
        // if excitation energy greater than ionisation threshold, then autoionisaiton
        if ((excitationEnergy > ioniThres) && (G4UniformRand() < 0.5)) {
          fParticleChangeForGamma->ProposeLocalEnergyDeposit(ioniThres);
          // energy of ejected electron
          G4double secondaryKinetic = excitationEnergy - ioniThres;
          // random direction
          G4double cosTheta = 2 * G4UniformRand() - 1., phi = CLHEP::twopi * G4UniformRand();
          G4double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
          G4double ux = sinTheta * std::cos(phi), uy = sinTheta * std::sin(phi), uz = cosTheta;
          G4ThreeVector deltaDirection(ux, uy, uz);
          // Create the new particle with its characteristics
          auto dp = new G4DynamicParticle(G4Electron::Electron(), deltaDirection, secondaryKinetic);
          fvect->push_back(dp);
        }
      }
      else {
        G4ExceptionDescription description;
        description << "Kinetic energy <= 0 at " << fpN2->GetName() << " material !!!";
        G4Exception("G4DNAPTBExcitationModel::SampleSecondaries", "", FatalException, description);
      }
    }
    else if (fpG4_WATER == nullptr || materialID != fpG4_WATER->GetIndex()) {
      // Retrieve the excitation energy for the current material
      G4double excitationEnergy = fTableMeanEnergyPTB[materialID];
      // Calculate the new energy of the particle
      G4double newEnergy = k - excitationEnergy;
      // Check that the new energy is above zero before applying it the particle.
      // Otherwise, do nothing.
      if (newEnergy > 0) {
        fParticleChangeForGamma->ProposeMomentumDirection(aDynamicParticle->GetMomentumDirection());
        fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
        fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
      }
      else {
        G4ExceptionDescription description;
        description << "Kinetic energy <= 0 at " << materialID << " index material !!!";
        G4Exception("G4DNAPTBExcitationModel::SampleSecondaries", "", FatalException, description);
      }
    }
    else {
      G4int level = RandomSelectShell(k, particle, materialID);
      G4double excitationEnergy = waterStructure.ExcitationEnergy(level);
      G4double newEnergy = k - excitationEnergy;

      if (newEnergy > 0) {
        fParticleChangeForGamma->ProposeMomentumDirection(aDynamicParticle->GetMomentumDirection());
        fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
        fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
        const G4Track* theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
        G4DNAChemistryManager::Instance()->CreateWaterMolecule(eExcitedMolecule, level,
                                                               theIncomingTrack);
      }
      else {
        G4ExceptionDescription description;
        description << "Kinetic energy <= 0 at " << materialID << " ID material !!!";
        G4Exception("G4DNAPTBExcitationModel::SampleSecondaries", "", FatalException, description);
      }
    }
  }
}
