
#include "Tst05PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

Tst05PrimaryGeneratorAction::Tst05PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
}

Tst05PrimaryGeneratorAction::~Tst05PrimaryGeneratorAction()
{
  delete particleGun;
}

void Tst05PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle 
    = particleTable->FindParticle(particleName="e-");    
  particleGun->SetParticleDefinition(particle);

  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(20.*MeV);  
  particleGun->SetParticlePosition(G4ThreeVector(0,0,0));

  particleGun->GeneratePrimaryVertex(anEvent);
}






