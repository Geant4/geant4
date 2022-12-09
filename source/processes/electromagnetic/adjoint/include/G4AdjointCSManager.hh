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
////////////////////////////////////////////////////////////////////////////////
//  Class:    G4AdjointCSManager
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
//
// Class is responsible for the management of all adjoint cross section
// matrices, and for the computation of the total forward and adjoint cross
// sections. Total adjoint and forward cross sections are needed to correct the
// weight of a particle after a tracking step or after the occurrence of a
// reverse reaction. It is also used to sample an adjoint secondary from a
// given adjoint cross section matrix.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef G4AdjointCSManager_h
#define G4AdjointCSManager_h 1

#include "globals.hh"
#include "G4AdjointCSMatrix.hh"
#include "G4ThreadLocalSingleton.hh"

#include <vector>

class G4Element;
class G4Material;
class G4MaterialCutsCouple;
class G4ParticleDefinition;
class G4PhysicsTable;
class G4VEmProcess;
class G4VEmAdjointModel;
class G4VEnergyLossProcess;

class G4AdjointCSManager
{
  friend class G4ThreadLocalSingleton<G4AdjointCSManager>;

 public:
  ~G4AdjointCSManager();
  static G4AdjointCSManager* GetAdjointCSManager();

  G4int GetNbProcesses();

  // Registration of the different models and processes

  std::size_t RegisterEmAdjointModel(G4VEmAdjointModel*);

  void RegisterEmProcess(G4VEmProcess* aProcess,
                         G4ParticleDefinition* aPartDef);

  void RegisterEnergyLossProcess(G4VEnergyLossProcess* aProcess,
                                 G4ParticleDefinition* aPartDef);

  void RegisterAdjointParticle(G4ParticleDefinition* aPartDef);

  // Building of the CS Matrices and Total Forward and Adjoint LambdaTables
  void BuildCrossSectionMatrices();

  void BuildTotalSigmaTables();

  // Get TotalCrossSections form Total Lambda Tables, Needed for Weight
  // correction and scaling of the
  G4double GetTotalAdjointCS(G4ParticleDefinition* aPartDef, G4double Ekin,
                             const G4MaterialCutsCouple* aCouple);

  G4double GetTotalForwardCS(G4ParticleDefinition* aPartDef, G4double Ekin,
                             const G4MaterialCutsCouple* aCouple);

  G4double GetAdjointSigma(G4double Ekin_nuc, std::size_t index_model,
                           G4bool is_scat_proj_to_proj,
                           const G4MaterialCutsCouple* aCouple);

  void GetEminForTotalCS(G4ParticleDefinition* aPartDef,
                         const G4MaterialCutsCouple* aCouple,
                         G4double& emin_adj, G4double& emin_fwd);

  void GetMaxFwdTotalCS(G4ParticleDefinition* aPartDef,
                        const G4MaterialCutsCouple* aCouple,
                        G4double& e_sigma_max, G4double& sigma_max);

  void GetMaxAdjTotalCS(G4ParticleDefinition* aPartDef,
                        const G4MaterialCutsCouple* aCouple,
                        G4double& e_sigma_max, G4double& sigma_max);

  // CrossSection Correction 1 or FwdCS/AdjCS following the G4boolean value of
  // forward_CS_is_used and forward_CS_mode
  G4double GetCrossSectionCorrection(G4ParticleDefinition* aPartDef,
                                     G4double PreStepEkin,
                                     const G4MaterialCutsCouple* aCouple,
                                     G4bool& fwd_is_used);

  // Cross section mode
  inline void SetFwdCrossSectionMode(G4bool aBool) { fForwardCSMode = aBool; }

  // Weight correction
  G4double GetContinuousWeightCorrection(G4ParticleDefinition* aPartDef,
                                         G4double PreStepEkin,
                                         G4double AfterStepEkin,
                                         const G4MaterialCutsCouple* aCouple,
                                         G4double step_length);

  G4double GetPostStepWeightCorrection();

  // called by the adjoint model to get the CS, if not otherwise specified
  G4double ComputeAdjointCS(G4Material* aMaterial, G4VEmAdjointModel* aModel,
                            G4double PrimEnergy, G4double Tcut,
                            G4bool isScatProjToProj,
                            std::vector<G4double>& AdjointCS_for_each_element);

  // called by the adjoint model to sample secondary energy from the CS matrix
  G4Element* SampleElementFromCSMatrices(G4Material* aMaterial,
                                         G4VEmAdjointModel* aModel,
                                         G4double PrimEnergy, G4double Tcut,
                                         G4bool isScatProjToProj);

  // Total Adjoint CS is computed at initialisation phase
  G4double ComputeTotalAdjointCS(const G4MaterialCutsCouple* aMatCutCouple,
                                 G4ParticleDefinition* aPart,
                                 G4double PrimEnergy);

  G4ParticleDefinition* GetAdjointParticleEquivalent(
    G4ParticleDefinition* theFwdPartDef);

  G4ParticleDefinition* GetForwardParticleEquivalent(
    G4ParticleDefinition* theAdjPartDef);

  // inline
  inline void SetIon(G4ParticleDefinition* adjIon, G4ParticleDefinition* fwdIon)
  {
    fAdjIon = adjIon;
    fFwdIon = fwdIon;
  }

 private:
  G4AdjointCSManager();

  void DefineCurrentMaterial(const G4MaterialCutsCouple* couple);

  void DefineCurrentParticle(const G4ParticleDefinition* aPartDef);

  G4double ComputeAdjointCS(G4double aPrimEnergy,
                            G4AdjointCSMatrix* anAdjointCSMatrix,
                            G4double Tcut);

  std::vector<G4AdjointCSMatrix*> BuildCrossSectionsModelAndElement(
    G4VEmAdjointModel* aModel, G4int Z, G4int A, G4int nbin_pro_decade);

  std::vector<G4AdjointCSMatrix*> BuildCrossSectionsModelAndMaterial(
    G4VEmAdjointModel* aModel, G4Material* aMaterial, G4int nbin_pro_decade);

  static constexpr G4double fTmin = 0.1 * CLHEP::keV;
  static constexpr G4double fTmax = 100. * CLHEP::TeV;
  // fNbins chosen to avoid error
  // in the CS value close to CS jump. (For example at Tcut)
  static constexpr G4int fNbins = 320;

  static G4ThreadLocal G4AdjointCSManager* fInstance;

  // only one ion can be considered by simulation
  G4ParticleDefinition* fAdjIon = nullptr;
  G4ParticleDefinition* fFwdIon = nullptr;

  G4MaterialCutsCouple* fCurrentCouple = nullptr;
  G4Material* fCurrentMaterial         = nullptr;

  // x dim is for G4VAdjointEM*, y dim is for elements
  std::vector<std::vector<G4AdjointCSMatrix*>>
    fAdjointCSMatricesForScatProjToProj;

  std::vector<std::vector<G4AdjointCSMatrix*>> fAdjointCSMatricesForProdToProj;

  std::vector<G4VEmAdjointModel*> fAdjointModels;

  std::vector<std::size_t> fIndexOfAdjointEMModelInAction;
  std::vector<G4bool> fIsScatProjToProj;
  std::vector<std::vector<G4double>> fLastAdjointCSVsModelsAndElements;

  // total adjoint and total forward cross section table in function of material
  // and in function of adjoint particle type
  std::vector<G4PhysicsTable*> fTotalFwdSigmaTable;
  std::vector<G4PhysicsTable*> fTotalAdjSigmaTable;

  // Sigma table for each G4VAdjointEMModel
  std::vector<G4PhysicsTable*> fSigmaTableForAdjointModelScatProjToProj;
  std::vector<G4PhysicsTable*> fSigmaTableForAdjointModelProdToProj;

  std::vector<std::vector<G4double>> fEminForFwdSigmaTables;
  std::vector<std::vector<G4double>> fEminForAdjSigmaTables;
  std::vector<std::vector<G4double>> fEkinofFwdSigmaMax;
  std::vector<std::vector<G4double>> fEkinofAdjSigmaMax;

  // list of forward G4VEmProcess and of G4VEnergyLossProcess for the different
  // adjoint particle
  std::vector<std::vector<G4VEmProcess*>*> fForwardProcesses;
  std::vector<std::vector<G4VEnergyLossProcess*>*> fForwardLossProcesses;

  // list of adjoint particles considered
  std::vector<G4ParticleDefinition*> fAdjointParticlesInAction;

  G4double fMassRatio              = 1.;  // ion
  G4double fLastCSCorrectionFactor = 1.;

  std::size_t fCurrentParticleIndex = 0;
  std::size_t fCurrentMatIndex      = 0;

  G4bool fCSMatricesBuilt = false;
  G4bool fSigmaTableBuilt = false;
  G4bool fForwardCSUsed   = true;
  G4bool fForwardCSMode   = true;
  // Two CS mode are possible:
  // 1) fForwardCSMode = false, the Adjoint CS are used as it is implying
  //    an AlongStep Weight Correction.
  // 2) fForwardCSMode = true, the Adjoint CS are scaled to have the total
  //    adjoint CS equal to the fwd one implying a PostStep Weight Correction.
  // For energies where the total Fwd CS or the total adjoint CS are zero,
  // the scaling is not possible and fForwardCSUsed is set to false
};
#endif
