#ifndef MySteppingAction_H
#define MySteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"


class MySteppingAction : public G4UserSteppingAction {

public:

  MySteppingAction();
  virtual ~MySteppingAction();

public:

  virtual void UserSteppingAction( const G4Step* );
  // The main method to define.

  inline G4double getTotalEdepAllParticles() const;
  // Total deposited energy, in the same event, 
  // due to all particles, in the whole experimental hall.

  inline void reset();
  // This method should be invoked at the end of the event,
  // to reset to zero the total deposit energy variables.

private:
 
  G4double totalEdepAllParticles;

};


inline G4double MySteppingAction::getTotalEdepAllParticles() const {
  return totalEdepAllParticles;
}


inline void MySteppingAction::reset() {
  totalEdepAllParticles = 0.0;
}


#endif


