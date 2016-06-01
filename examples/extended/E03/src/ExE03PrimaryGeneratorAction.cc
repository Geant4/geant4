
#include "ExE03PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

ExE03PrimaryGeneratorAction::ExE03PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
}

ExE03PrimaryGeneratorAction::~ExE03PrimaryGeneratorAction()
{
  delete particleGun;
}

void ExE03PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
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






