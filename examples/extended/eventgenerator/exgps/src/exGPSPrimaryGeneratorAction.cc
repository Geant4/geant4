
#include "exGPSPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

exGPSPrimaryGeneratorAction::exGPSPrimaryGeneratorAction()
{
   particleGun = new G4GeneralParticleSource();
}

exGPSPrimaryGeneratorAction::~exGPSPrimaryGeneratorAction()
{
  delete particleGun;
}

void exGPSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->GeneratePrimaryVertex(anEvent) ;
}






