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
// Class Description:
//
// This class is a modified clone of G4Cerenkov, extended to support
// offloading optical photon generation. Offloading can be enabled either
// via the G4OpticalParameters::Instance()->SetCerenkovOffloadPhotons(true)
// method or the UI command:
//
// /process/optical/cerenkov/setOffloadPhotons true
//
// When offloading is enabled, the process generates a single secondary track
// of type G4QuasiOpticalPhoton, along with associated metadata encapsulated
// in G4CerenkovQuasiTrackInfo. This auxiliary track information is used
// to generate optical photons at a later stage—typically during offloading.
//
// The intended workflow leverages G4VTrackingManager, which delegates these
// secondary tracks to a dedicated G4ProcessManager for G4QuasiOpticalPhoton.
// These tracks are then handled by a user-defined custom tracking manager,
// independent of the default process managers used for other particles.
//
// The primary purpose of this class is to facilitate the transfer of essential
// data for offloaded optical photon generation in heterogeneous computing
// models

#ifndef G4QuasiCerenkov_h
#define G4QuasiCerenkov_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4ForceCondition.hh"
#include "G4GPILSelection.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4VProcess.hh"

#include <map>

class G4Material;
class G4ParticleDefinition;
class G4PhysicsTable;
class G4Step;
class G4Track;
class G4VParticleChange;

class G4QuasiCerenkov : public G4VProcess
{
 public:
  explicit G4QuasiCerenkov(const G4String& processName = "QuasiCerenkov",
                      G4ProcessType type          = fElectromagnetic);
  ~G4QuasiCerenkov();

  explicit G4QuasiCerenkov(const G4QuasiCerenkov& right);

  G4QuasiCerenkov& operator=(const G4QuasiCerenkov& right) = delete;

  G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;
  // Returns true -> 'is applicable', for all charged particles
  // except short-lived particles.

  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) override;
  // Build table at a right time

  void PreparePhysicsTable(const G4ParticleDefinition& part) override;
  void Initialise();

  G4double GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*);
  // Returns the discrete step limit and sets the 'StronglyForced'
  // condition for the DoIt to be invoked at every step.

  G4double PostStepGetPhysicalInteractionLength(const G4Track& aTrack, G4double,
                                                G4ForceCondition*) override;
  // Returns the discrete step limit and sets the 'StronglyForced'
  // condition for the DoIt to be invoked at every step.

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                  const G4Step& aStep) override;
  // This is the method implementing the Cerenkov process.

  //  no operation in  AtRestDoIt and  AlongStepDoIt
  virtual G4double AlongStepGetPhysicalInteractionLength(
    const G4Track&, G4double, G4double, G4double&, G4GPILSelection*) override
  {
    return -1.0;
  };

  virtual G4double AtRestGetPhysicalInteractionLength(
    const G4Track&, G4ForceCondition*) override
  {
    return -1.0;
  };

  //  no operation in  AtRestDoIt and  AlongStepDoIt
  virtual G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) override
  {
    return nullptr;
  };

  virtual G4VParticleChange* AlongStepDoIt(const G4Track&,
                                           const G4Step&) override
  {
    return nullptr;
  };

  void SetTrackSecondariesFirst(const G4bool state);
  // If set, the primary particle tracking is interrupted and any
  // produced Cerenkov photons are tracked next. When all have
  // been tracked, the tracking of the primary resumes.

  G4bool GetTrackSecondariesFirst() const;
  // Returns the boolean flag for tracking secondaries first.

  void SetMaxBetaChangePerStep(const G4double d);
  // Set the maximum allowed change in beta = v/c in % (perCent) per step.

  G4double GetMaxBetaChangePerStep() const;
  // Returns the maximum allowed change in beta = v/c in % (perCent)

  void SetMaxNumPhotonsPerStep(const G4int NumPhotons);
  // Set the maximum number of Cerenkov photons allowed to be generated during
  // a tracking step. This is an average ONLY; the actual number will vary
  // around this average. If invoked, the maximum photon stack will roughly be
  // of the size set. If not called, the step is not limited by the number of
  // photons generated.

  G4int GetMaxNumPhotonsPerStep() const;
  // Returns the maximum number of Cerenkov photons allowed to be
  // generated during a tracking step.

  void SetStackPhotons(const G4bool);
  // Call by the user to set the flag for stacking the Cerenkov photons

  G4bool GetStackPhotons() const;
  // Return the boolean for whether or not the Cerenkov photons are stacked

  void SetOffloadPhotons(const G4bool);
  // Call by the user to set the flag for offloading the Cerenkov photons

  G4bool GetOffloadPhotons() const;
  // Return the boolean for whether or not the Cerenkov photons are offloaded

  G4int GetNumPhotons() const;
  // Returns the current number of scint. photons (after PostStepDoIt)

  G4PhysicsTable* GetPhysicsTable() const;
  // Returns the address of the physics table.

  void DumpPhysicsTable() const;
  // Prints the physics table.

  G4double GetAverageNumberOfPhotons(const G4double charge, const G4double beta,
                                     const G4Material* aMaterial,
                                     G4MaterialPropertyVector* Rindex) const;

  void DumpInfo() const override {ProcessDescription(G4cout);};
  void ProcessDescription(std::ostream& out) const override;

  void SetVerboseLevel(G4int);
  // sets verbosity

 protected:
  G4PhysicsTable* thePhysicsTable;
  std::map<std::size_t, std::size_t> fIndexMPT;

 private:
  G4double fMaxBetaChange;
  
  G4int fMaxPhotons;
  G4int fNumPhotons;

  G4bool fStackingFlag;
  G4bool fOffloadingFlag;
  G4bool fTrackSecondariesFirst;

  G4int secID = -1;  // creator modelID

};

inline G4bool G4QuasiCerenkov::GetTrackSecondariesFirst() const
{
  return fTrackSecondariesFirst;
}

inline G4double G4QuasiCerenkov::GetMaxBetaChangePerStep() const
{
  return fMaxBetaChange;
}

inline G4int G4QuasiCerenkov::GetMaxNumPhotonsPerStep() const { return fMaxPhotons; }

inline G4bool G4QuasiCerenkov::GetStackPhotons() const { return fStackingFlag; }

inline G4bool G4QuasiCerenkov::GetOffloadPhotons() const { return fOffloadingFlag; }

inline G4int G4QuasiCerenkov::GetNumPhotons() const { return fNumPhotons; }

inline G4PhysicsTable* G4QuasiCerenkov::GetPhysicsTable() const
{
  return thePhysicsTable;
}

#endif /* G4QuasiCerenkov_h */
