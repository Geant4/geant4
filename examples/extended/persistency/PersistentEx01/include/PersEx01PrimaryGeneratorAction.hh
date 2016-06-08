
#ifndef PersEx01PrimaryGeneratorAction_h
#define PersEx01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class PersEx01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PersEx01PrimaryGeneratorAction();
    ~PersEx01PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


