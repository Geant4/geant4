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
/////////////////////////////////////////////////////////////////////////////////
//  Class:    G4AdjointAlongStepWeightCorrection
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
/////////////////////////////////////////////////////////////////////////////////
//
//  Documentation:
//
//  Continuous processes act on adjoint particles to continuously correct their
//  weight during the adjoint reverse tracking. This process is needed when
//  the adjoint cross sections are not scaled such that the total adjoint cross
//  section matches the total forward cross section. By default the mode where
//  the total adjoint cross section is equal to the total forward cross section
//  is used and therefore this along step weightcorrection factor is 1. However
//  in some cases (some energy ranges) the total forward cross section or the
//  total adjoint cross section can be zero. In this case the along step weight
//  correction is needed and is given by exp(-(Sigma_tot_adj-Sigma_tot_fwd).dx)
//
//-------------------------------------------------------------

#ifndef G4AdjointAlongStepWeightCorrection_h
#define G4AdjointAlongStepWeightCorrection_h 1

#include "globals.hh"
#include "G4VContinuousProcess.hh"

class G4AdjointCSManager;
class G4MaterialCutsCouple;
class G4ParticleDefinition;
class G4ParticleChange;
class G4Step;
class G4Track;

class G4AdjointAlongStepWeightCorrection : public G4VContinuousProcess
{
 public:
  explicit G4AdjointAlongStepWeightCorrection(
    const G4String& name = "ContinuousWeightCorrection",
    G4ProcessType type   = fElectromagnetic);

  ~G4AdjointAlongStepWeightCorrection() override;

  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) override;

  void ProcessDescription(std::ostream&) const override;
  void DumpInfo() const override { ProcessDescription(G4cout); };

  G4AdjointAlongStepWeightCorrection(G4AdjointAlongStepWeightCorrection&) =
    delete;
  G4AdjointAlongStepWeightCorrection& operator=(
    const G4AdjointAlongStepWeightCorrection& right) = delete;

 protected:
  G4double GetContinuousStepLimit(const G4Track& track,
                                  G4double previousStepSize,
                                  G4double currentMinimumStep,
                                  G4double& currentSafety) override;

 private:
  void DefineMaterial(const G4MaterialCutsCouple* couple);

  const G4MaterialCutsCouple* fCurrentCouple = nullptr;
  G4AdjointCSManager* fCSManager = nullptr;
  G4ParticleChange* fParticleChange;

  G4double fPreStepKinEnergy = 1.;
};

inline void G4AdjointAlongStepWeightCorrection::DefineMaterial(
  const G4MaterialCutsCouple* couple)
{
  if(couple != fCurrentCouple)
  {
    fCurrentCouple = couple;
  }
}

#endif
