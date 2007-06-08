#ifndef DetectorSlicePrimaryGeneratorAction_h
#define DetectorSlicePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;


class DetectorSlicePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:

  DetectorSlicePrimaryGeneratorAction();
  ~DetectorSlicePrimaryGeneratorAction();

public:

  void GeneratePrimaries( G4Event* anEvent );

private:

    G4ParticleGun* particleGun;
};

#endif


