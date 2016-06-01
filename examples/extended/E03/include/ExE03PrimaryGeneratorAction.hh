
#ifndef ExE03PrimaryGeneratorAction_h
#define ExE03PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class ExE03PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ExE03PrimaryGeneratorAction();
    ~ExE03PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


