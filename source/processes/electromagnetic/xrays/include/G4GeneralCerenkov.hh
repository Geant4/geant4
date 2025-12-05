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
// --------------------------------------------------------------------
//
// G4GeneralCerenkov
//
// Class description:
// The Cerenkov process using model approach when an object G4VXRayModel
// is assign to G4LogicalVolume. A model is fully responsible for
// Cerenkov gamma production. This process class performing only tecnical
// operation and interaction with Geant4 kernel.
//
// Created 25.05.2025 V.Ivanchenko on base of G4Cerenkov class
//
// --------------------------------------------------------------------

#ifndef G4GeneralCerenkov_h
#define G4GeneralCerenkov_h 1

#include "globals.hh"
#include "G4ForceCondition.hh"
#include "G4LogicalVolume.hh"
#include "G4VXRayModel.hh"
#include "G4VDiscreteProcess.hh"

#include <vector>

class G4Material;
class G4ParticleDefinition;
class G4PhysicsTable;
class G4Step;
class G4Track;
class G4VParticleChange;

class G4GeneralCerenkov : public G4VDiscreteProcess
{
 public:
  explicit G4GeneralCerenkov(const G4String& processName = "Cerenkov",
                             G4ProcessType type = fElectromagnetic);
  ~G4GeneralCerenkov() override;

  G4GeneralCerenkov(const G4GeneralCerenkov& right) = delete;
  G4GeneralCerenkov& operator=(const G4GeneralCerenkov& right) = delete;

  G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;

  void PreparePhysicsTable(const G4ParticleDefinition& part) override;

  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) override;

  G4double PostStepGetPhysicalInteractionLength(const G4Track& aTrack, G4double,
                                                G4ForceCondition*) override;
  // Returns the discrete step limit and sets the 'StronglyForced'
  // condition for the DoIt to be invoked at every step.

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                  const G4Step& aStep) override;
  // This is the method implementing the Cerenkov process.

  void AddModelForVolume(G4VXRayModel*, const G4String& nameLogVolume);
  // explicit addition of a custom Cerenkov model to a logical volume

  void DumpInfo() const override { ProcessDescription(G4cout); };
  void ProcessDescription(std::ostream& out) const override;

  G4double GetMeanFreePath(const G4Track&, G4double,
			   G4ForceCondition*) override;
  
  // Obsolete methods to be removed for the next major release
  
  void SetTrackSecondariesFirst(const G4bool state);
  // If set, the primary particle tracking is interrupted and any
  // produced Cerenkov photons are tracked next. When all have
  // been tracked, the tracking of the primary resumes.

  void SetMaxBetaChangePerStep(const G4double d);
  // Set the maximum allowed change in beta = v/c in % (perCent) per step.

  void SetMaxNumPhotonsPerStep(const G4int NumPhotons);
  // Set the maximum number of Cerenkov photons allowed to be generated during
  // a tracking step. This is an average ONLY; the actual number will vary
  // around this average. If invoked, the maximum photon stack will roughly be
  // of the size set. If not called, the step is not limited by the number of
  // photons generated.

  void SetStackPhotons(const G4bool);
  // Call by the user to set the flag for stacking the Cerenkov photons

  void SetVerboseLevel(G4int);

private:

  const G4LogicalVolume* fCurrentLV{nullptr};
  G4VXRayModel* fCurrentModel{nullptr};

  G4double fMaxBetaChange{0.1};
  G4double fBetaMin{1.0};
  G4double fPreStepBeta{0.0};
  
  G4int fMaxPhotons{100};

  G4bool fStackingFlag{true};
  G4bool fTrackSecondariesFirst{true};
  G4bool isInitializer{false};
  G4bool isPrepared{false};
  G4bool isBuilt{false};

  G4int secID{-1};  // creator modelID
  G4int nModels{0};

  // map includes logical volume pointer and index of the model
  static std::vector<std::vector<const G4LogicalVolume*>* >* fLV;

  // vector is used only at initialisation
  // these models are destructed by G4LossTableManager
  static std::vector<G4VXRayModel*>* fSharedModels;

  // vector of names of logical volumes for master used for initilisation
  // not filled for a worker thread
  std::vector<G4String>* fLVNames{nullptr};
  
  // models used in run time - they are thread local, are
  // instantiated in worker thread, and are cloned from fSharedModels
  // these models are destructed by G4LossTableManager
  std::vector<G4VXRayModel*> fModels;

  // buffer for X-Rays
  std::vector<G4Track*> fSecondaries;
};

#endif
