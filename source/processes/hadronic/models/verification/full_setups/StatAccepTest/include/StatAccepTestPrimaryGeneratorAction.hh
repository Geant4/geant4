#ifndef StatAccepTestPrimaryGeneratorAction_h
#define StatAccepTestPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;


class StatAccepTestPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:

  StatAccepTestPrimaryGeneratorAction();
  ~StatAccepTestPrimaryGeneratorAction();

public:

  void GeneratePrimaries( G4Event* anEvent );

private:

    G4ParticleGun* particleGun;
};

#endif


