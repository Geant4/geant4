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
// G4DNADoubleIonisationModel.hh
//
//  Created at 2024/04/03 (Thu.)
//  Author: Shogo OKADA @KEK-CRC (shogo.okada@kek.jp)
//
//  Reference: J.Meesungnoen et. al, DOI: 10.1021/jp058037z
//

#ifndef G4DNA_DOUBLE_IONISATION_MODEL_HH_
#define G4DNA_DOUBLE_IONISATION_MODEL_HH_

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"

#include "G4DNAGenericIonsManager.hh"
#include "G4DNACrossSectionDataSet.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4LogLogInterpolation.hh"

#include "G4DNAWaterIonisationStructure.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4NistManager.hh"
#include "G4DNAMultipleIonisationManager.hh"

using EnergyLimitTable = std::map<G4String, G4double, std::less<G4String>>;

using CrossSectionDataTable = std::map<G4String, G4DNACrossSectionDataSet*,
                                       std::less<G4String>>;

//==============================================================================

class G4DNADoubleIonisationModel : public G4VEmModel {
public:

  // constructor
  G4DNADoubleIonisationModel(
    const G4ParticleDefinition* p = nullptr,
    const G4String& model_name = "G4DNADoubleIonisationModel");

  // destructor
  ~G4DNADoubleIonisationModel() override;

  G4DNADoubleIonisationModel& operator=(
                             const G4DNADoubleIonisationModel&) = delete;
  G4DNADoubleIonisationModel(const G4DNADoubleIonisationModel&) = delete;

  void Initialise(
    const G4ParticleDefinition* particle, const G4DataVector&) override;

  G4double CrossSectionPerVolume(
    const G4Material* material, const G4ParticleDefinition* pdef,
    G4double ekin, G4double, G4double) override;

  void SampleSecondaries(
    std::vector<G4DynamicParticle*>* vsec, const G4MaterialCutsCouple* couple,
    const G4DynamicParticle* particle, G4double, G4double) override;

  void SelectStationary(G4bool in);

  void SelectVerboseLevel(G4int in);

  void UseChampionAlphaParameter(G4bool in);

  void SetMultipleIonisationEnergy(G4double in);

protected:

  G4double RandomizeEjectedElectronEnergy(
    G4ParticleDefinition* pdef, G4double ekin, G4int shell);

  G4ParticleChangeForGamma* particle_change_ = nullptr;

  G4int RandomSelect(G4double energy, G4double scale_param,
                     const G4String& pname);

  G4double GenerateSecondaries(std::vector<G4DynamicParticle*>* vsec,
                               const G4MaterialCutsCouple* couple,
                               const G4DynamicParticle* particle,
                               G4int ioni_shell,
                               G4double& theta, G4double& phi,
                               G4double& shell_energy);

  G4double GetLowEnergyLimit(const G4String& pname);

  G4double GetUppEnergyLimit(const G4String& pname);

  G4bool stat_code_;

  G4VAtomDeexcitation* atom_deex_;

  EnergyLimitTable elow_tab_;
  EnergyLimitTable eupp_tab_;

  CrossSectionDataTable xs_tab_;

  G4ParticleDefinition* proton_def_{nullptr};
  G4ParticleDefinition* alpha_def_{nullptr};
  G4ParticleDefinition* carbon_def_{nullptr};

  const std::vector<G4double>* water_density_;

  G4bool is_initialized_;
  G4int verbose_level_;

  std::map<G4double, G4double> model_elow_tab_;

  G4DNAMultipleIonisationManager* mioni_manager_{nullptr};

  G4bool use_champion_param_;

  G4double energy_threshold_;
};

//==============================================================================
inline void G4DNADoubleIonisationModel::SelectStationary(G4bool in)
{
  stat_code_ = in;
}

//------------------------------------------------------------------------------
inline void G4DNADoubleIonisationModel::SelectVerboseLevel(G4int in)
{
  verbose_level_ = in;
}

//------------------------------------------------------------------------------
inline void G4DNADoubleIonisationModel::UseChampionAlphaParameter(G4bool in)
{
  use_champion_param_ = in;
}

//------------------------------------------------------------------------------
inline void G4DNADoubleIonisationModel::SetMultipleIonisationEnergy(G4double in)
{
  energy_threshold_ = in;
}

#endif // G4DNA_DOUBLE_IONISATION_MODEL_HH_
