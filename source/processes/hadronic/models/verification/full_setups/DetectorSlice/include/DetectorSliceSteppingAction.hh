#ifndef DetectorSliceSteppingAction_H
#define DetectorSliceSteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class DetectorSliceSteppingAction : public G4UserSteppingAction {

public:

  DetectorSliceSteppingAction();
  virtual ~DetectorSliceSteppingAction();

public:

  virtual void UserSteppingAction( const G4Step* );
  // The main method to define.

  inline G4double getTotalEdepAllParticlesInTracker() const;
  inline G4double getTotalEdepAllParticlesInEmCal() const;
  inline G4double getTotalEdepAllParticlesInHadCal() const;
  inline G4double getTotalEdepAllParticlesInMuonDetector() const;
  // Total deposited energy, in the same event, due to all particles,
  // inside, respectively, the following subdetectors: 
  // Tracker, Electromagnetic calorimeter, Hadronic calorimeter,
  // Muon detector.

  inline G4double getExitingRadiusPrimaryMuon() const;
  // Only in the case of a (primary) beam muon particle,
  // this method returns the muon exiting radius immediately
  // after the muon detector (in the case such primary muon
  // is still present).

  inline G4int getPrimaryParticleId() const;
  inline G4double getPrimaryParticleEnergy() const;
  // Id and Energy of the primary particle.

  void reset();
  // This method should be invoked at the end of the event,
  // to reset to zero the total deposit energy variables.

private:
 
  G4double totalEdepAllParticlesInTracker;
  G4double totalEdepAllParticlesInEmCal;
  G4double totalEdepAllParticlesInHadCal;
  G4double totalEdepAllParticlesInMuonDetector;

  G4double exitingRadiusPrimaryMuon;

  G4int    primaryParticleId;
  G4double primaryParticleEnergy;
  G4bool   isFirstStepOfTheEvent;

  G4bool isPrimaryMuonReachingMuonDetector;

  const G4int muonMinusId, muonPlusId;

};


inline G4double DetectorSliceSteppingAction::
getTotalEdepAllParticlesInTracker() const {
  return totalEdepAllParticlesInTracker;
}

inline G4double DetectorSliceSteppingAction::
getTotalEdepAllParticlesInEmCal() const {
  return totalEdepAllParticlesInEmCal;
}

inline G4double DetectorSliceSteppingAction::
getTotalEdepAllParticlesInHadCal() const {
  return totalEdepAllParticlesInHadCal;
}

inline G4double DetectorSliceSteppingAction::
getTotalEdepAllParticlesInMuonDetector() const {
  return totalEdepAllParticlesInMuonDetector;
}


inline G4double DetectorSliceSteppingAction::getExitingRadiusPrimaryMuon() const {
  return exitingRadiusPrimaryMuon;
}


inline G4int DetectorSliceSteppingAction::getPrimaryParticleId() const {
  return primaryParticleId;
}


inline G4double DetectorSliceSteppingAction::getPrimaryParticleEnergy() const {
  return primaryParticleEnergy;
}


#endif


