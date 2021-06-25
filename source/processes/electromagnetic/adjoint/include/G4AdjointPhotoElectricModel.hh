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
//  Class:   G4AdjointPhotoElectricModel
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
//
//  Model for the adjoint photo electric process.
//  Put a higher limit on the CS to avoid a high rate of Inverse Photo e-effect
//  at low energy. The very high adjoint CS of the reverse photo electric
//  reaction produce a high rate of reverse photo electric reaction in the inner
//  side of a shielding for eaxmple, the correction of this occurrence by weight
//  correction in the StepDoIt method is not statistically sufficient at small
//  energy. The problem is partially solved by setting a higher CS limit and
//  compensating it by an extra weight correction factor. However when coupling
//  it with other reverse processes the reverse photo-electric is still the
//  source of very occasional high weights that decrease the efficiency of the
//  computation. A way to solve this problemn is still needed but is difficult
//  to find as it happens in rare cases but does give a weight that is outside
//  the normal distribution. (Very Tricky!)
//
////////////////////////////////////////////////////////////////////////////////

#ifndef G4AdjointPhotoElectricModel_h
#define G4AdjointPhotoElectricModel_h 1

#include "globals.hh"
#include "G4VEmAdjointModel.hh"

class G4AdjointPhotoElectricModel : public G4VEmAdjointModel
{
 public:
  G4AdjointPhotoElectricModel();
  ~G4AdjointPhotoElectricModel() override;

  void SampleSecondaries(const G4Track& aTrack, G4bool isScatProjToProj,
                         G4ParticleChange* fParticleChange) override;

  G4double AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
                               G4double primEnergy,
                               G4bool isScatProjToProj) override;

  G4double AdjointCrossSectionPerAtom(const G4Element* anElement,
                                      G4double electronEnergy);

  G4AdjointPhotoElectricModel(G4AdjointPhotoElectricModel&) = delete;
  G4AdjointPhotoElectricModel& operator=(
    const G4AdjointPhotoElectricModel& right) = delete;

 protected:
  void CorrectPostStepWeight(G4ParticleChange* fParticleChange,
                             G4double old_weight, G4double adjointPrimKinEnergy,
                             G4double projectileKinEnergy,
                             G4bool isScatProjToProj) override;

 private:
  void DefineCurrentMaterialAndElectronEnergy(
    const G4MaterialCutsCouple* aCouple, G4double eEnergy);

  G4double fShellProb[40][40];
  G4double fXsec[40];
  G4double fTotAdjointCS      = 0.;
  G4double fFactorCSBiasing   = 1.;
  G4double fPreStepAdjointCS  = 0.;
  G4double fPostStepAdjointCS = 0.;
  G4double fCurrenteEnergy    = 0.;

  size_t fIndexElement = 0;
};

#endif
