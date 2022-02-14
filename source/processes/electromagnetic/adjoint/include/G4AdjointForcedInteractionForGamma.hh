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
//  Class:   G4AdjointForcedInteractionForGamma
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
//
//  Class for the forced interaction of reverse gamma
/////////////////////////////////////////////////////////////////////////////////

#ifndef G4AdjointForcedInteractionForGamma_h
#define G4AdjointForcedInteractionForGamma_h 1

#include "globals.hh"
#include "G4VContinuousDiscreteProcess.hh"

class G4VParticleChange;
class G4ParticleChange;
class G4ParticleDefinition;
class G4Track;
class G4VEmAdjointModel;
class G4AdjointCSManager;

class G4AdjointForcedInteractionForGamma : public G4VContinuousDiscreteProcess
{
 public:
  explicit G4AdjointForcedInteractionForGamma(G4String process_name);

  ~G4AdjointForcedInteractionForGamma() override;

  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  G4double PostStepGetPhysicalInteractionLength(
    const G4Track& track, G4double previousStepSize,
    G4ForceCondition* condition) override;

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

  G4VParticleChange* AlongStepDoIt(const G4Track& track,
                                   const G4Step& step) override;

  inline void RegisterAdjointComptonModel(G4VEmAdjointModel* adjModel)
  {
    fAdjointComptonModel = adjModel;
  }

  inline void RegisterAdjointBremModel(G4VEmAdjointModel* adjModel)
  {
    fAdjointBremModel = adjModel;
  }

  void ProcessDescription(std::ostream&) const override;
  void DumpInfo() const override { ProcessDescription(G4cout); };

  G4AdjointForcedInteractionForGamma(G4AdjointForcedInteractionForGamma&) =
    delete;
  G4AdjointForcedInteractionForGamma& operator=(
    const G4AdjointForcedInteractionForGamma& right) = delete;

 protected:
  G4double GetMeanFreePath(const G4Track& track, G4double previousStepSize,
                           G4ForceCondition* condition) override;

  G4double GetContinuousStepLimit(const G4Track& aTrack,
                                  G4double previousStepSize,
                                  G4double currentMinimumStep,
                                  G4double& currentSafety) override;

 private:
  G4VEmAdjointModel* fAdjointComptonModel;
  G4VEmAdjointModel* fAdjointBremModel;

  G4ParticleChange* fParticleChange;
  G4AdjointCSManager* fCSManager;

  G4double fLastAdjCS = 0.;
  G4double fCSBias = 1.;
  G4double fAccTrackLength    = 0.;
  G4double fTotNbAdjIntLength = 0.;

  G4double fNbAdjIntLength = 0.;

  G4bool fContinueGammaAsNewFreeFlight = false;
  G4bool fFreeFlightGamma              = false;
  G4bool fCopyGammaForForced           = false;
};

#endif
