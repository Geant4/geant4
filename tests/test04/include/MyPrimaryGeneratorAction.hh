
#ifndef MyPrimaryGeneratorAction_h
#define MyPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class MyPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    MyPrimaryGeneratorAction();
    ~MyPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


