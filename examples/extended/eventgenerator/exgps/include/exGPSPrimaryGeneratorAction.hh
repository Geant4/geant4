
#ifndef exGPSPrimaryGeneratorAction_h
#define exGPSPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4GeneralParticleSource;
class G4Event;

class exGPSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    exGPSPrimaryGeneratorAction();
    ~exGPSPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:

    G4GeneralParticleSource* particleGun;
};

#endif


