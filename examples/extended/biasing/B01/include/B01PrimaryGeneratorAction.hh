
#ifndef B01PrimaryGeneratorAction_hh
#define B01PrimaryGeneratorAction_hh B01PrimaryGeneratorAction_hh 

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class B01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    B01PrimaryGeneratorAction();
    ~B01PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


