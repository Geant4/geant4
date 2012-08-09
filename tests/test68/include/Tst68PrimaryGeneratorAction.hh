#ifndef Tst68PrimaryGeneratorAction_h
#define Tst68PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;


class Tst68PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:

  Tst68PrimaryGeneratorAction();
  ~Tst68PrimaryGeneratorAction();

public:

  void GeneratePrimaries( G4Event* anEvent );

private:

    G4ParticleGun* particleGun;
};

#endif
