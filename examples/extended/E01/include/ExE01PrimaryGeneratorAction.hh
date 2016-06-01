
#ifndef ExE01PrimaryGeneratorAction_h
#define ExE01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class ExE01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ExE01PrimaryGeneratorAction();
    ~ExE01PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


