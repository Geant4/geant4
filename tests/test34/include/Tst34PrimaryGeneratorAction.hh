#ifndef Tst34PrimaryGeneratorAction_h
#define Tst34PrimaryGeneratorAction_h

#include "G4ThreeVector.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4Event;

class Tst34PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:
	Tst34PrimaryGeneratorAction();
	~Tst34PrimaryGeneratorAction();
	void GeneratePrimaries(G4Event* anEvent);

private:
  G4GeneralParticleSource	* particleGun;
};

#endif


