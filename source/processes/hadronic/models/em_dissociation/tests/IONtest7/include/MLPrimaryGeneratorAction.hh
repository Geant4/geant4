#ifndef MLPrimaryGeneratorAction_h
#define MLPrimaryGeneratorAction_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"
////////////////////////////////////////////////////////////////////////////////
//
class MLPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  MLPrimaryGeneratorAction ();
  ~MLPrimaryGeneratorAction ();

public:
  void GeneratePrimaries (G4Event* anEvent);
  G4GeneralParticleSource* GetParticleGun () const {return particleGun;};

private:
  G4GeneralParticleSource* particleGun;
  
};
////////////////////////////////////////////////////////////////////////////////
#endif
