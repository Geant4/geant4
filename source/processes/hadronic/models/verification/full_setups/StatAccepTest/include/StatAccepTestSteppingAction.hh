#ifndef StatAccepTestSteppingAction_H
#define StatAccepTestSteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class StatAccepTestSteppingAction : public G4UserSteppingAction {

public:

  StatAccepTestSteppingAction();
  virtual ~StatAccepTestSteppingAction();

public:

  virtual void UserSteppingAction( const G4Step* );
  // The main method to define.

  inline G4double getTotalEdepAllParticles() const;
  // Total deposited energy, in the same event, due to all particles,
  // inside the calorimeter.

  inline G4int getPrimaryParticleId() const;
  inline G4double getPrimaryParticleEnergy() const;
  // Id and Energy of the primary particle.

  void reset();
  // This method should be invoked at the end of the event,
  // to reset to zero the total deposit energy variables.

private:
 
  G4double totalEdepAllParticles;
  G4int    primaryParticleId;
  G4double primaryParticleEnergy;
  G4bool   isFirstStepOfTheEvent;

};


inline G4double StatAccepTestSteppingAction::getTotalEdepAllParticles() const {
  return totalEdepAllParticles;
}


inline G4int StatAccepTestSteppingAction::getPrimaryParticleId() const {
  return primaryParticleId;
}


inline G4double StatAccepTestSteppingAction::getPrimaryParticleEnergy() const {
  return primaryParticleEnergy;
}


#endif


