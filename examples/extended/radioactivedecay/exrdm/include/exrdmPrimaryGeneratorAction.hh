

//
// **********************************************************************

#ifndef exrdmPrimaryGeneratorAction_h
#define exrdmPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4GeneralParticleSource;
class G4Event;

class exrdmPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  exrdmPrimaryGeneratorAction();
  ~exrdmPrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);

private:
  G4GeneralParticleSource* particleGun;
  
};

#endif



