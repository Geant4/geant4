
#ifndef Tst05PrimaryGeneratorAction_h
#define Tst05PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class Tst05PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Tst05PrimaryGeneratorAction();
    ~Tst05PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


