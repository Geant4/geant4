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
// G4DNADoubleIonisationModel.cc
//
//  Created at 2024/04/03 (Thu.)
//  Author: Shogo OKADA @KEK-CRC (shogo.okada@kek.jp)
//
//  Reference: J.Meesungnoen et. al, DOI: 10.1021/jp058037z
//

#include "G4DNADoubleIonisationModel.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

#include "G4IonTable.hh"
#include "G4GenericIon.hh"
#include "G4DNARuddAngle.hh"
#include "G4DeltaAngle.hh"
#include "G4Exp.hh"

#include <sstream>

namespace {

G4DNAWaterIonisationStructure water_structure;

// parameters for rejection function
struct FuncParams {
  G4double Bj_energy;
  G4double alpha_const;
  G4double beta_squared;
  G4double velocity;
  G4double correction_factor;
  G4double wc;
  G4double F1;
  G4double F2;
  G4double c;
};

//------------------------------------------------------------------------------
void setup_rejection_function(G4ParticleDefinition* pdef, const G4double ekin,
                              const G4int shell, FuncParams& par)
{

  // Following values provided by M. Dingfelder (priv. comm)
  const G4double Bj[5]
    = { 12.60 * eV, 14.70 * eV, 18.40 * eV, 32.20 * eV, 540.0 * eV };

  // Data For Liquid Water from Dingfelder (Protons in Water)
  G4double A1{1.02}, B1{82.0}, C1{0.45}, D1{-0.80}, E1{0.38}, A2{1.07},
           B2{11.6}, // Value provided by M. Dingfelder (priv. comm)
           C2{0.60}, D2{0.04}, alpha_const{0.64};

  auto Bj_energy = Bj[shell];

  if (shell == 4) {
    alpha_const = 0.66;
    //Data For Liquid Water K SHELL from Dingfelder (Protons in Water)
    A1 = 1.25; B1 = 0.5;  C1 = 1.00; D1 = 1.00; E1 = 3.00;
    A2 = 1.10; B2 = 1.30; C2 = 1.00; D2 = 0.00;
    // The following cases are provided by M. Dingfelder (priv. comm)
    Bj_energy = water_structure.IonisationEnergy(shell);
  }

  const auto mass  = pdef->GetPDGMass();
  const auto tau   = ekin * electron_mass_c2 / mass;
  const auto A_ion = pdef->GetAtomicMass();

  G4double v2;
  G4double beta2;

  constexpr G4double Ry  = 13.6 * eV;
  constexpr G4double xxx = 5.447761194E-02 * MeV;

  if (tau < xxx) {
    v2 = tau / Bj_energy;
    beta2 = 2.0 * tau / electron_mass_c2;
  } else {
    // Relativistic
    v2 = (0.5 * electron_mass_c2 / Bj_energy)
            * (1.0 - (1.0 / std::pow((1.0 + (tau / electron_mass_c2)), 2.0)));
    beta2 = 1.0 - 1.0 / std::pow((1.0 + (tau / electron_mass_c2 / A_ion)), 2.0);
  }

  const auto v = std::sqrt(v2);
  const auto wc = 4.0 * v2 - 2.0 * v - (Ry / (4.0 * Bj_energy));

  const auto L1 = (C1 * std::pow(v, D1)) / (1.0 + E1 * std::pow(v, (D1 + 4.0)));
  const auto L2 = C2 * std::pow(v, D2);
  const auto H1 = (A1 * G4Log(1.0 + v2)) / (v2 + (B1 / v2));
  const auto H2 = (A2 / v2) + (B2 /(v2 * v2));
  const auto F1 = L1 + H1;
  const auto F2 = (L2 * H2) / (L2 + H2);

  // ZF. generalized & relativistic version
  G4double max_energy;

  if (ekin <= 0.1 * mass) {
    // maximum kinetic energy , non relativistic
    max_energy = 4.0 * (electron_mass_c2 / mass) * ekin;
  } else {
    // relativistic
    auto gamma = 1.0 / std::sqrt(1.0 - beta2);
    max_energy = 2.0 * electron_mass_c2 * (gamma * gamma - 1.0)
        / (1.0 + 2.0 * gamma * (electron_mass_c2 / mass)
            + std::pow(electron_mass_c2 / mass, 2.0));
  }

  const auto wmax = max_energy / Bj_energy;
  auto c = wmax * (F2 * wmax+ F1 * (2.0 + wmax))
              / (2.0 * (1.0 + wmax) * (1.0 + wmax));
  c = 1.0 / c; // manual calculus leads to c = 1 / c

  par.Bj_energy = Bj_energy;
  par.alpha_const = alpha_const;
  par.beta_squared = beta2;
  par.velocity = v;
  par.correction_factor = 1.0;
  par.wc = wc;
  par.F1 = F1;
  par.F2 = F2;
  par.c = c;

}

//------------------------------------------------------------------------------
G4double rejection_function(G4ParticleDefinition* pdef, const G4int shell,
                            const FuncParams& par, G4double proposed_ws)
{

  const G4double Gj[5] = { 0.99, 1.11, 1.11, 0.52, 1.0 };

  proposed_ws /= par.Bj_energy;

  auto rejection_term = 1.0 + G4Exp(par.alpha_const * (proposed_ws - par.wc)
                                      / par.velocity);
  rejection_term = (1.0 / rejection_term) * par.correction_factor * Gj[shell];

  if (pdef == G4Proton::ProtonDefinition()) {

    // for protons
    return rejection_term;

  } else if (pdef->GetAtomicMass() > 4) {

    // for carbon ions
    auto Z = pdef->GetAtomicNumber();
    auto x = 100.0 * std::sqrt(par.beta_squared) / std::pow(Z, 0.6666667);
    auto zeff = Z * (1.0 - G4Exp(x * (-1.316 + x * (0.112 - 0.0650 * x))));

    rejection_term *= (zeff * zeff);

    return rejection_term;

  }

  // for alpha particles
  auto zeff = pdef->GetPDGCharge() / eplus + pdef->GetLeptonNumber();
  rejection_term *= (zeff * zeff);

  return rejection_term;

}

//------------------------------------------------------------------------------
G4double proposed_sampled_energy(const FuncParams& par)
{
  const auto rval = G4UniformRand();

  auto proposed_ws = par.c * (par.F1 * par.F1 * par.c
                      + 2.0 * rval * (par.F2 - par.F1));

  proposed_ws = -par.F1 * par.c + 2.0 * rval + std::sqrt(proposed_ws);

  proposed_ws /= (par.c * (par.F1 + par.F2) - 2.0 * rval);

  proposed_ws *= par.Bj_energy;

  return proposed_ws;
}

} // end of anonymous namespace

//==============================================================================

// constructor
G4DNADoubleIonisationModel::G4DNADoubleIonisationModel(
  const G4ParticleDefinition*, const G4String& model_name)
    : G4VEmModel(model_name),
      is_initialized_(false)
{
  water_density_ = nullptr;

  model_elow_tab_[1] = 100 * eV;
  model_elow_tab_[4] = 1.0 * keV;
  model_elow_tab_[5] = 0.5 * MeV; // For A = 3 or above, limit is MeV/uma

  verbose_level_ = 0;

  // Define default angular generator
  SetAngularDistribution(new G4DNARuddAngle());

  // Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);
  atom_deex_ = nullptr;
  particle_change_ = nullptr;

  // Selection of stationary mode
  stat_code_ = false;

  // True if use champion alpha parameter
  use_champion_param_ = false;

  // Double-ionization energy
  energy_threshold_ = 40.0 * eV;
}

//------------------------------------------------------------------------------
G4DNADoubleIonisationModel::~G4DNADoubleIonisationModel()
{
  for (auto x : xs_tab_) {
    G4DNACrossSectionDataSet* table = x.second;
    if (table) { delete table; }
  }
}

//------------------------------------------------------------------------------
void G4DNADoubleIonisationModel::Initialise(
  const G4ParticleDefinition* particle, const G4DataVector&)
{
  if (verbose_level_ > 3) {
    G4cout << "Calling G4DNADoubleIonisationModel::Initialise()" << G4endl;
  }

  proton_def_  = G4Proton::ProtonDefinition();
  alpha_def_   = G4DNAGenericIonsManager::Instance()->GetIon("alpha++");
  carbon_def_  = G4IonTable::GetIonTable()->GetIon(6, 12);

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
    G4cout << "G4DNADoubleIonisationModel is initialized " << G4endl
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
G4double G4DNADoubleIonisationModel::GetLowEnergyLimit(const G4String& pname)
{
  G4double elim{0.0};

  EnergyLimitTable::iterator itr = elow_tab_.find(pname);
  if (itr != elow_tab_.end()) { elim = itr->second; }

  return elim;
}

//------------------------------------------------------------------------------
G4double G4DNADoubleIonisationModel::GetUppEnergyLimit(const G4String& pname)
{
  G4double elim{0.0};

  EnergyLimitTable::iterator itr = eupp_tab_.find(pname);
  if (itr != eupp_tab_.end()) { elim = itr->second; }

  return elim;
}

//------------------------------------------------------------------------------
G4double G4DNADoubleIonisationModel::CrossSectionPerVolume(
  const G4Material* material, const G4ParticleDefinition* pdef,
  G4double ekin, G4double, G4double)
{

  if (verbose_level_ > 3) {
    G4cout << "Calling G4DNADoubleIonisationModel::CrossSectionPerVolume()"
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
      G4Exception("G4DNADoubleIonisationModel::CrossSectionPerVolume",
                  "em0002", FatalException,
                  "Model not applicable to particle type.");
    }

    G4DNACrossSectionDataSet* table = pos->second;
    if (table != nullptr) {
      const auto a = mioni_manager_->GetAlphaParam(ekin);
      sigma = table->FindValue(ekin) * a;
    }

  }

  if (verbose_level_ > 2) {

    std::stringstream msg;

    msg << "----------------------------------------------------------------\n";
    msg << " G4DNADoubleIonisationModel - XS INFO START\n";
    msg << "  - Kinetic energy(eV): " << ekin/eV << ", Particle : "
        << pdef->GetParticleName() << "\n";
    msg << "  - Cross section per water molecule (cm^2):  "
        << sigma / cm / cm << "\n";
    msg << "  - Cross section per water molecule (cm^-1): "
        << sigma * water_dens / (1.0 / cm) << "\n";
    msg << " G4DNADoubleIonisationModel - XS INFO END\n";
    msg << "----------------------------------------------------------------\n";

    G4cout << msg.str() << G4endl;

  }

  return (sigma * water_dens);

}

//------------------------------------------------------------------------------
G4double G4DNADoubleIonisationModel::GenerateSecondaries(
  std::vector<G4DynamicParticle*>* vsec, const G4MaterialCutsCouple* couple,
  const G4DynamicParticle* particle, G4int ioni_shell,
  G4double& theta, G4double& phi, G4double& shell_energy)
{
  auto pdef = particle->GetDefinition();

  // get kinetic energy for a parent particle
  auto ekin1 = particle->GetKineticEnergy();

  // sample kinetic energy for a secondary electron
  auto ekin2 = RandomizeEjectedElectronEnergy(pdef, ekin1, ioni_shell);

  // sample momentum direction for a secondary electron
  auto sample_electron_direction = [this](
      const G4DynamicParticle* dp, G4double _ekin2, G4int _Z, G4int _ioni_shell,
      const G4MaterialCutsCouple* mcc, G4double& _theta, G4double& _phi) {

    G4ThreeVector locdir;

    if (_theta > 0.0) {

      auto costh = std::cos(_theta);
      auto sinth = std::sqrt((1.0 - costh) * (1.0 + costh));
      locdir.set(sinth * std::cos(_phi), sinth * std::sin(_phi), costh);
      locdir.rotateUz(dp->GetMomentumDirection());

    } else {

      locdir = GetAngularDistribution()->SampleDirectionForShell(
                  dp, _ekin2, _Z, _ioni_shell, mcc->GetMaterial());
      _theta = locdir.theta();
      _phi   = locdir.phi();

    }

    return locdir;
  };

  constexpr G4int Z = 8;
  auto delta_dir = sample_electron_direction(
                     particle, ekin2, Z, ioni_shell, couple, theta, phi);

  // generate a secondary electron and put it into the stack
  auto dp = new G4DynamicParticle(G4Electron::Electron(), delta_dir, ekin2);
  vsec->push_back(dp);

  if (!atom_deex_ || ioni_shell != 4) { return ekin2; }

  // ***************************************************************************
  // Only atomic deexcitation from K shell is considered

  constexpr auto k_shell = G4AtomicShellEnumerator(0);
  const auto shell = atom_deex_->GetAtomicShell(Z, k_shell);

  // get number of secondary electrons in the stack
  // before processing atomic deescitation
  const auto num_sec_init = vsec->size();

  // perform atomic deexcitation process
  atom_deex_->GenerateParticles(vsec, shell, Z, 0, 0);

  // get number of secondary electrons in the stack
  // after processing atomic deescitation
  const auto num_sec_final = vsec->size();

  if (num_sec_final == num_sec_init) { return ekin2; }

  for (auto i = num_sec_init; i < num_sec_final; i++) {

    auto e = ((*vsec)[i])->GetKineticEnergy();

    // Check if there is enough residual energy
    if (shell_energy < e) {

      // Invalid secondary: not enough energy to create it!
      // Keep its energy in the local deposit
      delete (*vsec)[i];
      (*vsec)[i] = 0;

      continue;

    }

    // Ok, this is a valid secondary: keep it
    shell_energy -= e;
  }

  // ***************************************************************************

  return ekin2;
}

//------------------------------------------------------------------------------
void G4DNADoubleIonisationModel::SampleSecondaries(
  std::vector<G4DynamicParticle*>* vsec, const G4MaterialCutsCouple* couple,
  const G4DynamicParticle* particle, G4double, G4double)
{

  if (verbose_level_ > 3) {
    G4cout << "Calling SampleSecondaries() of G4DNADoubleIonisationModel"
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

  constexpr G4int kNumSecondaries = 2;
  constexpr G4double kDeltaTheta  = pi;

  G4int ioni_shell[kNumSecondaries];
  G4double shell_energy[kNumSecondaries];

  const auto scale_param = mioni_manager_->GetAlphaParam(ekin);
  G4double tot_ioni_energy{0.0};
  for (G4int i = 0; i < kNumSecondaries; i++) {
    ioni_shell[i]   = RandomSelect(ekin, scale_param, pname);
    shell_energy[i] = ::water_structure.IonisationEnergy(ioni_shell[i]);
    tot_ioni_energy += shell_energy[i];
  }

  if (ekin < tot_ioni_energy || tot_ioni_energy < energy_threshold_) {
    return;
  }

  // generate secondary electrons
  G4double theta{0.0}, phi{0.0}, tot_ekin2{0.0};
  for (G4int i = 0; i < kNumSecondaries; i++) {
    tot_ekin2 += GenerateSecondaries(vsec, couple, particle, ioni_shell[i],
                                     theta, phi, shell_energy[i]);
    theta += kDeltaTheta;
  }

  // This should never happen
  if (mioni_manager_->CheckShellEnergy(eDoubleIonisedMolecule, shell_energy)) {
    G4Exception("G4DNADoubleIonisatioModel::SampleSecondaries()",
        "em2050", FatalException, "Negative local energy deposit");
  }

  // ***************************************************************************
  // update kinematics for this parent particle
  const auto primary_dir = particle->GetMomentumDirection();
  particle_change_->ProposeMomentumDirection(primary_dir);

  const auto scattered_energy = ekin - tot_ioni_energy - tot_ekin2;

  // update total amount of shell energy
  tot_ioni_energy = shell_energy[0] + shell_energy[1];

  if (stat_code_) {
    particle_change_->SetProposedKineticEnergy(ekin);
    particle_change_->ProposeLocalEnergyDeposit(ekin - scattered_energy);
  } else {
    particle_change_->SetProposedKineticEnergy(scattered_energy);
    particle_change_->ProposeLocalEnergyDeposit(tot_ioni_energy);
  }

  // ***************************************************************************
  // generate double-ionized water molecules (H2O^2+)
  const auto the_track = particle_change_->GetCurrentTrack();
  mioni_manager_->CreateMultipleIonisedWaterMolecule(
                    eDoubleIonisedMolecule, ioni_shell, the_track);
  // ***************************************************************************

}

//------------------------------------------------------------------------------
G4double G4DNADoubleIonisationModel::RandomizeEjectedElectronEnergy(
  G4ParticleDefinition* pdef, G4double ekin, G4int shell)
{

  //
  // based on RandomizeEjectedElectronEnergy()
  //     of G4DNARuddIonisationExtendedModel
  //

  ::FuncParams par;
  ::setup_rejection_function(pdef, ekin, shell, par);

  // calculate maximum value
  G4double emax{0.0}, val;
  for (G4double en = 0.0; en < 20.0; en += 1.0) {
    val = ::rejection_function(pdef, shell, par, en);
    if (val <= emax) { continue; }
    emax = val;
  }

  G4double proposed_energy, rand;
  do {
    // Proposed energy by inverse function sampling
    proposed_energy = ::proposed_sampled_energy(par);
    rand = G4UniformRand() * emax;
    val = ::rejection_function(pdef, shell, par, proposed_energy);
  } while (rand > val);

  return proposed_energy;
}

//------------------------------------------------------------------------------
G4int G4DNADoubleIonisationModel::RandomSelect(
  G4double ekin, G4double scale_param, const G4String& pname)
{

  //
  // based on RandomSelect() of G4DNARuddIonisationExtendedModel
  //

  // Retrieve data table corresponding to the current particle type
  CrossSectionDataTable::iterator pos = xs_tab_.find(pname);

  if (pos == xs_tab_.end()) {
    G4Exception("G4DNADoubleIonisationModel::RandomSelect", "em0002",
        FatalException, "Model not applicable to particle type.");
  }

  G4DNACrossSectionDataSet* table = pos->second;

  if (table != nullptr) {

    // get total number of energy level
    const auto num_component = table->NumberOfComponents();

    auto* valuesBuffer = new G4double[num_component];

    auto shell = num_component;
    G4double value = 0.0;

    while (shell > 0) {
      shell--;
      valuesBuffer[shell] = table->GetComponent((G4int)shell)->FindValue(ekin)
                              * scale_param;
      value += valuesBuffer[shell];
    }

    value *= G4UniformRand();

    shell = num_component;

    while (shell > 0) {
      shell--;
      if (valuesBuffer[shell] > value) {
        if (valuesBuffer) { delete [] valuesBuffer; }
        return (G4int)shell;
      }
      value -= valuesBuffer[shell];
    }

    if (valuesBuffer) { delete [] valuesBuffer; }

  }

  return 0;
}
