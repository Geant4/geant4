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
// G4DNATripleIonisationModel.cc
//
//  Created at 2024/04/03 (Thu.)
//  Author: Shogo OKADA @KEK-CRC (shogo.okada@kek.jp)
//
//  Reference: J.Meesungnoen et. al, DOI: 10.1021/jp058037z
//

#include "G4DNATripleIonisationModel.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4IonTable.hh"
#include "G4GenericIon.hh"
#include "G4DNARuddAngle.hh"

#include <sstream>

namespace {

G4DNAWaterIonisationStructure water_structure;

} // end of anonymous namespace

//==============================================================================

// constructor
G4DNATripleIonisationModel::G4DNATripleIonisationModel(
  const G4ParticleDefinition* p, const G4String& model_name)
    : G4DNADoubleIonisationModel(p, model_name)
{
  // Triple-ionisation energy
  energy_threshold_ = 65.0 * eV;
}

//------------------------------------------------------------------------------
void G4DNATripleIonisationModel::Initialise(
  const G4ParticleDefinition* particle, const G4DataVector&)
{
  if (verbose_level_ > 3) {
    G4cout << "Calling G4DNATripleIonisationModel::Initialise()" << G4endl;
  }

  proton_def_ = G4Proton::ProtonDefinition();
  alpha_def_  = G4DNAGenericIonsManager::Instance()->GetIon("alpha++");
  carbon_def_ = G4IonTable::GetIonTable()->GetIon(6, 12);

  constexpr G4double kScaleFactor = 1.0 * m * m;

  mioni_manager_ = new G4DNAMultipleIonisationManager();

  G4double Z{0.0}, A{0.0};
  G4String alpha_param_file{"dna/multipleionisation_alphaparam_champion.dat"};

  if (particle == proton_def_) {

    // *************************************************************************
    // for protons
    auto proton = proton_def_->GetParticleName();
    elow_tab_[proton] = model_elow_tab_[1];
    eupp_tab_[proton] = 3.0 * MeV;

    // load cross-section data for single ionization process
    auto xs_proton = new G4DNACrossSectionDataSet(
                          new G4LogLogInterpolation, eV, kScaleFactor);
    xs_proton->LoadData("dna/sigma_ionisation_p_rudd");
    xs_tab_[proton] = xs_proton;

    // set energy limits
    SetLowEnergyLimit(elow_tab_[proton]);
    SetHighEnergyLimit(eupp_tab_[proton]);

    if (!use_champion_param_) {
      alpha_param_file = "dna/multipleionisation_alphaparam_p.dat";
    }

    Z = static_cast<G4double>(proton_def_->GetAtomicNumber());
    A = static_cast<G4double>(proton_def_->GetAtomicMass());

  } else if (particle == alpha_def_) {

    //**************************************************************************
    // for alpha particles
    auto alpha = alpha_def_->GetParticleName();
    elow_tab_[alpha] = model_elow_tab_[4];
    eupp_tab_[alpha] = 23.0 * MeV;

    // load cross-section data for single ionization process
    auto xs_alpha = new G4DNACrossSectionDataSet(
                        new G4LogLogInterpolation, eV, kScaleFactor);
    xs_alpha->LoadData("dna/sigma_ionisation_alphaplusplus_rudd");
    xs_tab_[alpha] = xs_alpha;

    // set energy limits
    SetLowEnergyLimit(elow_tab_[alpha]);
    SetHighEnergyLimit(eupp_tab_[alpha]);

    if (!use_champion_param_) {
      alpha_param_file = "dna/multipleionisation_alphaparam_alphaplusplus.dat";
    }

    Z = static_cast<G4double>(alpha_def_->GetAtomicNumber());
    A = static_cast<G4double>(alpha_def_->GetAtomicMass());

  } else if (particle == G4GenericIon::GenericIonDefinition()) {

    // *************************************************************************
    // for carbon ions
    auto carbon = carbon_def_->GetParticleName();
    elow_tab_[carbon] = model_elow_tab_[5] * carbon_def_->GetAtomicMass();
    eupp_tab_[carbon] = 120.0 * MeV;

    // load cross-section data for single ionization process
    auto xs_carbon = new G4DNACrossSectionDataSet(
                        new G4LogLogInterpolation, eV, kScaleFactor);
    xs_carbon->LoadData("dna/sigma_ionisation_c_rudd");
    xs_tab_[carbon] = xs_carbon;

    // set energy limits
    SetLowEnergyLimit(elow_tab_[carbon]);
    SetHighEnergyLimit(eupp_tab_[carbon]);

    if (!use_champion_param_) {
      alpha_param_file = "dna/multipleionisation_alphaparam_c.dat";
    }

    Z = static_cast<G4double>(carbon_def_->GetAtomicNumber());
    A = static_cast<G4double>(carbon_def_->GetAtomicMass());

  }

  // load alpha parameter
  mioni_manager_->LoadAlphaParam(alpha_param_file, Z, A);

  if (verbose_level_ > 0) {
    G4cout << "G4DNATripleIonisationModel is initialized " << G4endl
           << "Energy range: "
           << LowEnergyLimit() / eV << " eV - "
           << HighEnergyLimit() / keV << " keV for "
           << particle->GetParticleName()
           << G4endl;
  }

  water_density_ = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(
                      G4Material::GetMaterial("G4_WATER"));

  atom_deex_ = G4LossTableManager::Instance()->AtomDeexcitation();

  if (is_initialized_) { return; }

  particle_change_ = GetParticleChangeForGamma();
  is_initialized_ = true;
}

//------------------------------------------------------------------------------
G4double G4DNATripleIonisationModel::CrossSectionPerVolume(
  const G4Material* material, const G4ParticleDefinition* pdef,
  G4double ekin, G4double, G4double)
{

  if (verbose_level_ > 3) {
    G4cout << "Calling G4DNATripleIonisationModel::CrossSectionPerVolume()"
           << G4endl;
  }

  // Calculate total cross section for model

  if (pdef != proton_def_ && pdef != alpha_def_ && pdef != carbon_def_) {
    return 0.0;
  }

  static G4double water_dens = (*water_density_)[material->GetIndex()];

  const auto& pname = pdef->GetParticleName();

  const auto low_energy_lim = GetLowEnergyLimit(pname);
  const auto upp_energy_lim = GetUppEnergyLimit(pname);

  G4double sigma{0.0};
  if (ekin <= upp_energy_lim) {

    if (ekin < low_energy_lim) { ekin = low_energy_lim; }

    CrossSectionDataTable::iterator pos = xs_tab_.find(pname);
    if (pos == xs_tab_.end()) {
      G4Exception("G4DNATripleIonisationModel::CrossSectionPerVolume",
                  "em0002", FatalException,
                  "Model not applicable to particle type.");
    }

    G4DNACrossSectionDataSet* table = pos->second;
    if (table != nullptr) {
      const auto a = mioni_manager_->GetAlphaParam(ekin);
      sigma = table->FindValue(ekin) * a * a;
    }

  }

  if (verbose_level_ > 2) {

    std::stringstream msg;

    msg << "----------------------------------------------------------------\n";
    msg << " G4DNATripleIonisationModel - XS INFO START\n";
    msg << "  - Kinetic energy(eV): " << ekin/eV << ", Particle : "
        << pdef->GetParticleName() << "\n";
    msg << "  - Cross section per water molecule (cm^2):  "
        << sigma / cm / cm << "\n";
    msg << "  - Cross section per water molecule (cm^-1): "
        << sigma * water_dens / (1.0 / cm) << "\n";
    msg << " G4DNATripleIonisationModel - XS INFO END\n";
    msg << "----------------------------------------------------------------\n";

    G4cout << msg.str() << G4endl;

  }

  return (sigma * water_dens);

}

//------------------------------------------------------------------------------
void G4DNATripleIonisationModel::SampleSecondaries(
  std::vector<G4DynamicParticle*>* vsec, const G4MaterialCutsCouple* couple,
  const G4DynamicParticle* particle, G4double, G4double)
{

  if (verbose_level_ > 3) {
    G4cout << "Calling SampleSecondaries() of G4DNATripleIonisationModel"
           << G4endl;
  }

  // get the definition for this parent particle
  auto pdef = particle->GetDefinition();

  // get kinetic energy
  auto ekin = particle->GetKineticEnergy();

  // get particle name
  const auto& pname = pdef->GetParticleName();

  // get energy limits
  const auto low_energy_lim = GetLowEnergyLimit(pname);

  // ***************************************************************************
  // stop the transportation process of this parent particle
  // if its kinetic energy  is below the lower limit
  if (ekin < low_energy_lim) {
    particle_change_->SetProposedKineticEnergy(0.0);
    particle_change_->ProposeTrackStatus(fStopAndKill);
    particle_change_->ProposeLocalEnergyDeposit(ekin);
    return;
  }
  // ***************************************************************************

  constexpr G4int kNumSecondaries = 3;
  constexpr G4double kDeltaTheta  = pi * 0.666666667;

  G4int ioni_shell[kNumSecondaries] = {0, 0, 0};
  G4double shell_energy[kNumSecondaries];

  auto scale_param = mioni_manager_->GetAlphaParam(ekin);
  scale_param *= scale_param;

  G4bool is_continue{true};
  while (1) {
    ioni_shell[0] = RandomSelect(ekin, scale_param, pname);
    ioni_shell[1] = RandomSelect(ekin, scale_param, pname);
    ioni_shell[2] = RandomSelect(ekin, scale_param, pname);
    is_continue = (ioni_shell[0] == ioni_shell[1] &&
                   ioni_shell[1] == ioni_shell[2]);
    if (!is_continue) { break; }
  }

  G4double tot_ioni_energy{0.0};
  for (int i = 0; i < kNumSecondaries; i++) {
    shell_energy[i] = ::water_structure.IonisationEnergy(ioni_shell[i]);
    tot_ioni_energy += shell_energy[i];
  }

  if (ekin < tot_ioni_energy || tot_ioni_energy < energy_threshold_) {
    return;
  }

  // generate secondary electrons
  G4double theta{0.0}, phi{0.0}, tot_ekin2{0.0};
  for (int i = 0; i < kNumSecondaries; i++) {
    tot_ekin2 += GenerateSecondaries(vsec, couple, particle, ioni_shell[i],
                                     theta, phi, shell_energy[i]);
    theta += kDeltaTheta;
  }

  // This should never happen
  if (mioni_manager_->CheckShellEnergy(eTripleIonisedMolecule, shell_energy)) {
    G4Exception("G4DNATripleIonisatioModel::SampleSecondaries()",
        "em2050", FatalException, "Negative local energy deposit");
  }

  // ***************************************************************************
  // update kinematics for this parent particle
  const auto primary_dir = particle->GetMomentumDirection();
  particle_change_->ProposeMomentumDirection(primary_dir);

  const auto scattered_energy = ekin - tot_ioni_energy - tot_ekin2;

  // update total amount of shell energy
  tot_ioni_energy = shell_energy[0] + shell_energy[1] + shell_energy[2];

  if (stat_code_) {
    particle_change_->SetProposedKineticEnergy(ekin);
    particle_change_->ProposeLocalEnergyDeposit(
                                ekin - scattered_energy);
  } else {
    particle_change_->SetProposedKineticEnergy(scattered_energy);
    particle_change_->ProposeLocalEnergyDeposit(tot_ioni_energy);
  }

  // ***************************************************************************
  // generate triple-ionized water molecules (H2O^3+)
  const auto the_track = particle_change_->GetCurrentTrack();
  mioni_manager_->CreateMultipleIonisedWaterMolecule(
                    eTripleIonisedMolecule, ioni_shell, the_track);
  // ***************************************************************************

}
