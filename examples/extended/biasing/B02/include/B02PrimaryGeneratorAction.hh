
#ifndef B02PrimaryGeneratorAction_hh
#define B02PrimaryGeneratorAction_hh B02PrimaryGeneratorAction_hh 

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class B02PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    B02PrimaryGeneratorAction();
    ~B02PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


