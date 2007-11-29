#include "ExP02PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

ExP02PrimaryGeneratorAction::ExP02PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  //    particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="pi+"));
  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="geantino"));
  particleGun->SetParticleEnergy(1.0*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, 0.0));
}

ExP02PrimaryGeneratorAction::~ExP02PrimaryGeneratorAction()
{
  delete particleGun;
}

void ExP02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ThreeVector v(1.0,0.1,0.1);

  particleGun->SetParticleMomentumDirection(v);
  particleGun->GeneratePrimaryVertex(anEvent);
}


