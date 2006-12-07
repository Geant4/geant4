#ifndef RandomCaloSteppingAction_H
#define RandomCaloSteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"


class RandomCaloSteppingAction : public G4UserSteppingAction {

public:

  RandomCaloSteppingAction();
  virtual ~RandomCaloSteppingAction();

public:

  virtual void UserSteppingAction( const G4Step* );
  // The main method to define.

  inline G4double getTotalEdepAllParticles() const;
  // Total deposited energy, in the same event, due to all particles,
  // inside the calorimeter.

  void reset();
  // This method should be invoked at the end of the event,
  // to reset to zero the total deposit energy variables.

private:
 
  G4double totalEdepAllParticles;
  G4bool   isFirstStepOfTheEvent;

};


inline G4double RandomCaloSteppingAction::getTotalEdepAllParticles() const {
  return totalEdepAllParticles;
}


#endif


