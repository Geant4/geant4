
#ifndef LXePrimaryGeneratorAction_h
#define LXePrimaryGeneratorAction_h 1

#include "G4GeneralParticleSource.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class LXePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  LXePrimaryGeneratorAction();
  ~LXePrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);
  
private:
    G4ParticleGun* particleGun;
};

#endif


