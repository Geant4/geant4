
#ifndef ExE02PrimaryGeneratorAction_h
#define ExE02PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class ExE02PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ExE02PrimaryGeneratorAction();
    ~ExE02PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


