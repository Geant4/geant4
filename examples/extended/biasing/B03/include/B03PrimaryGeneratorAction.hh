
#ifndef B03PrimaryGeneratorAction_hh
#define B03PrimaryGeneratorAction_hh B03PrimaryGeneratorAction_hh 

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class B03PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    B03PrimaryGeneratorAction();
    ~B03PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


