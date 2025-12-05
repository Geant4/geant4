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
// This class is a modified clone of G4Scintillation, extended to support
// offloading optical photon generation. Offloading can be enabled either
// via the G4OpticalParameters::Instance()->SetScintOffloadPhotons(true)
// method or the UI command:
//
// /process/optical/scintillation/setOffloadPhotons true
//
// When offloading is enabled, the process generates a single secondary track
// of type G4QuasiOpticalPhoton, along with associated metadata encapsulated
// in G4ScintillationQuasiTrackInfo. This auxiliary track information is used
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

#ifndef G4QuasiScintillation_h
#define G4QuasiScintillation_h 1

#include "globals.hh"
#include "G4EmSaturation.hh"
#include "G4OpticalPhoton.hh"
#include "G4VRestDiscreteProcess.hh"

#include <map>

class G4PhysicsTable;
class G4Step;
class G4Track;

class G4QuasiScintillation : public G4VRestDiscreteProcess
{
 public:
  explicit G4QuasiScintillation(const G4String& procName = "QausiScintillation",
                                G4ProcessType type       = fElectromagnetic);
  ~G4QuasiScintillation();

  G4QuasiScintillation(const G4QuasiScintillation& right) = delete;
  G4QuasiScintillation& operator=(const G4QuasiScintillation& right) = delete;

  // G4QuasiScintillation Process has both PostStepDoIt (for energy
  // deposition of particles in flight) and AtRestDoIt (for energy
  // given to the medium by particles at rest)

  G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;
  // Returns true -> 'is applicable', for any particle type except
  // for an 'opticalphoton' and for short-lived particles

  void ProcessDescription(std::ostream&) const override;
  void DumpInfo() const override {ProcessDescription(G4cout);};

  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) override;
  // Build table at the right time

  void PreparePhysicsTable(const G4ParticleDefinition& part) override;
  void Initialise();

  G4double GetMeanFreePath(const G4Track& aTrack, G4double,
                           G4ForceCondition*) override;
  // Returns infinity; i. e. the process does not limit the step,
  // but sets the 'StronglyForced' condition for the DoIt to be
  // invoked at every step.

  G4double GetMeanLifeTime(const G4Track& aTrack, G4ForceCondition*) override;
  // Returns infinity; i. e. the process does not limit the time,
  // but sets the 'StronglyForced' condition for the DoIt to be
  // invoked at every step.

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                  const G4Step& aStep) override;
  G4VParticleChange* AtRestDoIt(const G4Track& aTrack,
                                const G4Step& aStep) override;

  G4double GetScintillationYieldByParticleType(const G4Track& aTrack,
                                               const G4Step& aStep,
                                               G4double& yield1,
                                               G4double& yield2,
                                               G4double& yield3,
                                               G4double& timeconstant1,
                                               G4double& timeconstant2,
                                               G4double& timeconstant3);
  // allow multiple time constants with scint by particle type
  // Returns the number of scintillation photons calculated when
  // scintillation depends on the particle type and energy
  // deposited (includes nonlinear dependendency) and updates the
  // yields for each channel

  void SetTrackSecondariesFirst(const G4bool state);
  // If set, the primary particle tracking is interrupted and any
  // produced scintillation photons are tracked next. When all
  // have been tracked, the tracking of the primary resumes.

  G4bool GetTrackSecondariesFirst() const;
  // Returns the boolean flag for tracking secondaries first.

  void SetFiniteRiseTime(const G4bool state);
  // If set, the G4QuasiScintillation process expects the user to have
  // set the constant material property SCINTILLATIONRISETIME{1,2,3}.

  G4bool GetFiniteRiseTime() const;
  // Returns the boolean flag for a finite scintillation rise time.

  G4PhysicsTable* GetIntegralTable1() const;
  // Returns the address of scintillation integral table #1.

  G4PhysicsTable* GetIntegralTable2() const;
  // Returns the address of scintillation integral table #2.

  G4PhysicsTable* GetIntegralTable3() const;
  // Returns the address of scintillation integral table #3.

  void AddSaturation(G4EmSaturation* sat);
  // Adds Birks Saturation to the process.

  void RemoveSaturation();
  // Removes the Birks Saturation from the process.

  G4EmSaturation* GetSaturation() const;
  // Returns the Birks Saturation.

  void SetScintillationByParticleType(const G4bool);
  // Called by the user to set the scintillation yield as a function
  // of energy deposited by particle type

  G4bool GetScintillationByParticleType() const;
  // Return the boolean that determines the method of scintillation
  // production

  void SetScintillationTrackInfo(const G4bool trackType);
  // Call by the user to set the G4ScintillationTrackInformation
  // to scintillation photon track

  G4bool GetScintillationTrackInfo() const;
  // Return the boolean for whether or not the
  // G4QuasiScintillationTrackInformation is set to the scint. photon track

  void SetStackPhotons(const G4bool);
  // Call by the user to set the flag for stacking the scint. photons

  G4bool GetStackPhotons() const;
  // Return the boolean for whether or not the scint. photons are stacked

  void SetOffloadPhotons(const G4bool);
  // Call by the user to set the flag for offloading the scint. photons

  G4bool GetOffloadPhotons() const;
  // Return the boolean for whether or not the scint. photons are offloaded

  G4int GetNumPhotons() const;
  // Returns the current number of scint. photons (after PostStepDoIt)

  void DumpPhysicsTable() const;
  // Prints the fast and slow scintillation integral tables.

  void SetVerboseLevel(G4int);
  // sets verbosity

 private:
  void BuildInverseCdfTable(const G4MaterialPropertyVector* MPV,
                            G4PhysicsFreeVector* vec) const;
  // Build the inverse cumulative distribution function (C.D.F.) table
  // for the scintillation photon energy spectrum

 private:

  G4PhysicsTable* fIntegralTable1;
  G4PhysicsTable* fIntegralTable2;
  G4PhysicsTable* fIntegralTable3;
  std::map<std::size_t, std::size_t> fIndexMPT;

  G4EmSaturation* fEmSaturation;
  const G4ParticleDefinition* opticalphoton =
    G4OpticalPhoton::OpticalPhotonDefinition();

  G4int fNumPhotons;
  
  G4bool fScintillationByParticleType;
  G4bool fScintillationTrackInfo;
  G4bool fStackingFlag;
  G4bool fOffloadingFlag;
  G4bool fTrackSecondariesFirst;
  G4bool fFiniteRiseTime;

#ifdef G4DEBUG_SCINTILLATION
  G4double ScintTrackEDep, ScintTrackYield;
#endif
  // emission time distribution when there is a finite rise time
  G4double sample_time(G4double tau1, G4double tau2);

  G4int secID = -1;  // creator modelID
  G4int fNumEnergyWarnings = 0;

};

////////////////////
// Inline methods
////////////////////

inline G4bool G4QuasiScintillation::GetTrackSecondariesFirst() const
{
  return fTrackSecondariesFirst;
}

inline G4bool G4QuasiScintillation::GetFiniteRiseTime() const
{
  return fFiniteRiseTime;
}

inline G4PhysicsTable* G4QuasiScintillation::GetIntegralTable1() const
{
  return fIntegralTable1;
}

inline G4PhysicsTable* G4QuasiScintillation::GetIntegralTable2() const
{
  return fIntegralTable2;
}

inline G4PhysicsTable* G4QuasiScintillation::GetIntegralTable3() const
{
  return fIntegralTable3;
}

inline void G4QuasiScintillation::AddSaturation(G4EmSaturation* sat)
{
  fEmSaturation = sat;
}

inline void G4QuasiScintillation::RemoveSaturation()
{
  fEmSaturation = nullptr;
}

inline G4EmSaturation* G4QuasiScintillation::GetSaturation() const
{
  return fEmSaturation;
}

inline G4bool G4QuasiScintillation::GetScintillationByParticleType() const
{
  return fScintillationByParticleType;
}

inline G4bool G4QuasiScintillation::GetScintillationTrackInfo() const
{
  return fScintillationTrackInfo;
}

inline G4bool G4QuasiScintillation::GetStackPhotons() const
{
  return fStackingFlag;
}

inline G4bool G4QuasiScintillation::GetOffloadPhotons() const
{
  return fOffloadingFlag;
}

inline G4int G4QuasiScintillation::GetNumPhotons() const
{
  return fNumPhotons;
}

#endif /* G4QuasiScintillation_h */
