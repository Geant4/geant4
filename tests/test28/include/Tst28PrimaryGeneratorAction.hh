#ifndef Tst28PrimaryGeneratorAction_h
#define Tst28PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class Tst28PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Tst28PrimaryGeneratorAction();
    ~Tst28PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


