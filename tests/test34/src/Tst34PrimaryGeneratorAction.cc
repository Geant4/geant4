#include "Tst34PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4ios.hh"

Tst34PrimaryGeneratorAction::Tst34PrimaryGeneratorAction()
{
	particleGun=new G4GeneralParticleSource;
}

Tst34PrimaryGeneratorAction::~Tst34PrimaryGeneratorAction()
{
	delete particleGun;
}

void Tst34PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
	particleGun->GeneratePrimaryVertex(anEvent);
}
