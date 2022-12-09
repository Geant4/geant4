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
//  Class:  G4VEMAdjointModel
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
//
//  Base class for Adjoint EM model. It is based on the use of direct
//  G4VEmModel.
////////////////////////////////////////////////////////////////////////////////

#ifndef G4VEmAdjointModel_h
#define G4VEmAdjointModel_h 1

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmModel.hh"

class G4AdjointCSMatrix;
class G4AdjointCSManager;
class G4Material;
class G4MaterialCutsCouple;
class G4ParticleChange;
class G4Region;
class G4Track;

class G4VEmAdjointModel
{
 public:
  explicit G4VEmAdjointModel(const G4String& nam);

  virtual ~G4VEmAdjointModel();

  //------------------------------------------------------------------------
  // Virtual methods to be implemented for the sample secondaries concrete model
  //------------------------------------------------------------------------

  virtual void SampleSecondaries(const G4Track& aTrack, G4bool isScatProjToProj,
                                 G4ParticleChange* fParticleChange) = 0;

  //------------------------------------------------------------------------
  // Methods for adjoint processes
  //------------------------------------------------------------------------

  virtual G4double AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
                                       G4double primEnergy,
                                       G4bool isScatProjToProj);

  // The implementation of the DiffCrossSection... here are correct for
  // energy loss process. For the photoelectric and Compton scattering
  // the method should be redefined
  virtual G4double DiffCrossSectionPerAtomPrimToSecond(
    G4double kinEnergyProj,  // kin energy of primary before interaction
    G4double kinEnergyProd,  // kinetic energy of the secondary particle
    G4double Z, G4double A = 0.);

  virtual G4double DiffCrossSectionPerAtomPrimToScatPrim(
    G4double kinEnergyProj,      // kin energy of primary before interaction
    G4double kinEnergyScatProj,  // kin energy of primary after interaction
    G4double Z, G4double A = 0.);

  virtual G4double DiffCrossSectionPerVolumePrimToSecond(
    const G4Material* aMaterial,
    G4double kinEnergyProj,  // kin energy of primary before interaction
    G4double kinEnergyProd   // kinetic energy of secondary particle
  );

  virtual G4double DiffCrossSectionPerVolumePrimToScatPrim(
    const G4Material* aMaterial,
    G4double kinEnergyProj,     // kin energy of primary before interaction
    G4double kinEnergyScatProj  // kinetic energy of primary after interaction
  );

  // Energy limits of adjoint secondary
  //------------------

  virtual G4double GetSecondAdjEnergyMaxForScatProjToProj(
    G4double primAdjEnergy);

  virtual G4double GetSecondAdjEnergyMinForScatProjToProj(
    G4double primAdjEnergy, G4double tcut = 0.);

  virtual G4double GetSecondAdjEnergyMaxForProdToProj(G4double primAdjEnergy);

  virtual G4double GetSecondAdjEnergyMinForProdToProj(G4double primAdjEnergy);

  // Other Methods
  //---------------

  void DefineCurrentMaterial(const G4MaterialCutsCouple* couple);

  std::vector<std::vector<double>*>
  ComputeAdjointCrossSectionVectorPerAtomForSecond(G4double kinEnergyProd,
                                                   G4double Z, G4double A = 0.,
                                                   G4int nbin_pro_decade = 10);

  std::vector<std::vector<double>*>
  ComputeAdjointCrossSectionVectorPerAtomForScatProj(
    G4double kinEnergyProd, G4double Z, G4double A = 0.,
    G4int nbin_pro_decade = 10);

  std::vector<std::vector<double>*>
  ComputeAdjointCrossSectionVectorPerVolumeForSecond(
    G4Material* aMaterial, G4double kinEnergyProd, G4int nbin_pro_decade = 10);

  std::vector<std::vector<double>*>
  ComputeAdjointCrossSectionVectorPerVolumeForScatProj(
    G4Material* aMaterial, G4double kinEnergyProd, G4int nbin_pro_decade = 10);

  inline void SetCSMatrices(std::vector<G4AdjointCSMatrix*>* Vec1CSMatrix,
                            std::vector<G4AdjointCSMatrix*>* Vec2CSMatrix)
  {
    fCSMatrixProdToProjBackScat = Vec1CSMatrix;
    fCSMatrixProjToProjBackScat = Vec2CSMatrix;
  };

  inline G4ParticleDefinition*
  GetAdjointEquivalentOfDirectPrimaryParticleDefinition()
  {
    return fAdjEquivDirectPrimPart;
  }

  inline G4ParticleDefinition*
  GetAdjointEquivalentOfDirectSecondaryParticleDefinition()
  {
    return fAdjEquivDirectSecondPart;
  }

  inline G4double GetHighEnergyLimit() { return fHighEnergyLimit; }

  inline G4double GetLowEnergyLimit() { return fLowEnergyLimit; }

  void SetHighEnergyLimit(G4double aVal);

  void SetLowEnergyLimit(G4double aVal);

  inline void DefineDirectEMModel(G4VEmModel* aModel) { fDirectModel = aModel; }

  void SetAdjointEquivalentOfDirectPrimaryParticleDefinition(
    G4ParticleDefinition* aPart);

  inline void SetAdjointEquivalentOfDirectSecondaryParticleDefinition(
    G4ParticleDefinition* aPart)
  {
    fAdjEquivDirectSecondPart = aPart;
  }

  inline void SetSecondPartOfSameType(G4bool aBool)
  {
    fSecondPartSameType = aBool;
  }

  inline G4bool GetSecondPartOfSameType() { return fSecondPartSameType; }

  inline void SetUseMatrix(G4bool aBool) { fUseMatrix = aBool; }

  inline void SetUseMatrixPerElement(G4bool aBool)
  {
    fUseMatrixPerElement = aBool;
  }

  inline void SetUseOnlyOneMatrixForAllElements(G4bool aBool)
  {
    fOneMatrixForAllElements = aBool;
  }

  inline void SetApplyCutInRange(G4bool aBool) { fApplyCutInRange = aBool; }

  inline G4bool GetUseMatrix() { return fUseMatrix; }

  inline G4bool GetUseMatrixPerElement() { return fUseMatrixPerElement; }

  inline G4bool GetUseOnlyOneMatrixForAllElements()
  {
    return fOneMatrixForAllElements;
  }

  inline G4bool GetApplyCutInRange() { return fApplyCutInRange; }

  inline G4String GetName() { return fName; }

  inline virtual void SetCSBiasingFactor(G4double aVal)
  {
    fCsBiasingFactor = aVal;
  }

  inline void SetCorrectWeightForPostStepInModel(G4bool aBool)
  {
    fInModelWeightCorr = aBool;
  }

  inline void SetAdditionalWeightCorrectionFactorForPostStepOutsideModel(
    G4double factor)
  {
    fOutsideWeightFactor = factor;
  }

  G4VEmAdjointModel(G4VEmAdjointModel&) = delete;
  G4VEmAdjointModel& operator=(const G4VEmAdjointModel& right) = delete;

 protected:
  G4double DiffCrossSectionFunction1(G4double kinEnergyProj);

  G4double DiffCrossSectionFunction2(G4double kinEnergyProj);

  // General methods to sample secondary energy
  G4double SampleAdjSecEnergyFromCSMatrix(std::size_t MatrixIndex,
                                          G4double prim_energy,
                                          G4bool isScatProjToProj);

  G4double SampleAdjSecEnergyFromCSMatrix(G4double prim_energy,
                                          G4bool isScatProjToProj);

  void SelectCSMatrix(G4bool isScatProjToProj);

  virtual G4double SampleAdjSecEnergyFromDiffCrossSectionPerAtom(
    G4double prim_energy, G4bool isScatProjToProj);

  // Post  Step weight correction
  virtual void CorrectPostStepWeight(G4ParticleChange* fParticleChange,
                                     G4double old_weight,
                                     G4double adjointPrimKinEnergy,
                                     G4double projectileKinEnergy,
                                     G4bool isScatProjToProj);

  G4AdjointCSManager* fCSManager;
  G4VEmModel* fDirectModel = nullptr;

  const G4String fName;

  G4Material* fSelectedMaterial        = nullptr;
  G4Material* fCurrentMaterial         = nullptr;
  G4MaterialCutsCouple* fCurrentCouple = nullptr;

  // particle definition
  G4ParticleDefinition* fAdjEquivDirectPrimPart   = nullptr;
  G4ParticleDefinition* fAdjEquivDirectSecondPart = nullptr;
  G4ParticleDefinition* fDirectPrimaryPart        = nullptr;

  // adjoint CS matrix for each element or material
  std::vector<G4AdjointCSMatrix*>* fCSMatrixProdToProjBackScat = nullptr;
  std::vector<G4AdjointCSMatrix*>* fCSMatrixProjToProjBackScat = nullptr;

  std::vector<G4double> fElementCSScatProjToProj;
  std::vector<G4double> fElementCSProdToProj;

  G4double fKinEnergyProdForIntegration     = 0.;
  G4double fKinEnergyScatProjForIntegration = 0.;

  G4double fLastCS                         = 0.;
  G4double fLastAdjointCSForScatProjToProj = 0.;
  G4double fLastAdjointCSForProdToProj     = 0.;

  G4double fPreStepEnergy = 0.;

  G4double fTcutPrim   = 0.;
  G4double fTcutSecond = 0.;

  // Energy limits
  G4double fHighEnergyLimit = 0.;
  G4double fLowEnergyLimit  = 0.;

  // Cross Section biasing factor
  G4double fCsBiasingFactor = 1.;

  // [1] This is needed for the forced interaction where part of the weight
  // correction is given outside the model while the secondary are created in
  // the model. The weight should be fixed before adding the secondary
  G4double fOutsideWeightFactor = 1.;

  // Needed for CS integration at the initialisation phase
  G4int fASelectedNucleus = 0;
  G4int fZSelectedNucleus = 0;

  std::size_t fCSMatrixUsed = 0;  // Index of crosssection matrices used

  G4bool fSecondPartSameType = false;
  G4bool fInModelWeightCorr =
    false;  // correct_weight_for_post_step_in_model, see [1]

  G4bool fApplyCutInRange = true;

  // Type of Model with Matrix or not
  G4bool fUseMatrix               = false;
  G4bool fUseMatrixPerElement     = false;  // other possibility is per Material
  G4bool fOneMatrixForAllElements = false;
};

#endif
