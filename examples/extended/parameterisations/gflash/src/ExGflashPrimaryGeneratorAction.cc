#include "ExGflashPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4ios.hh"

ExGflashPrimaryGeneratorAction::ExGflashPrimaryGeneratorAction()
{
	particleGun=new G4GeneralParticleSource;
}

ExGflashPrimaryGeneratorAction::~ExGflashPrimaryGeneratorAction()
{
	delete particleGun;
}

void ExGflashPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
	particleGun->GeneratePrimaryVertex(anEvent);
}
