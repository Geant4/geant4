#ifndef ExGflashPrimaryGeneratorAction_h
#define ExGflashPrimaryGeneratorAction_h

#include "G4ThreeVector.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4Event;

class ExGflashPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:
	ExGflashPrimaryGeneratorAction();
	~ExGflashPrimaryGeneratorAction();
	void GeneratePrimaries(G4Event* anEvent);

private:
  G4GeneralParticleSource	* particleGun;
};

#endif


