// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "ExN06PrimaryGeneratorAction.hh"

#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

ExN06PrimaryGeneratorAction::ExN06PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

}

ExN06PrimaryGeneratorAction::~ExN06PrimaryGeneratorAction()
{
  delete particleGun;
}

void ExN06PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="e+");

  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleTime(0.0*ns);

  G4ThreeVector position(0.0*cm,0.0*cm,0.0*cm);
  particleGun->SetParticlePosition(position);

  G4ThreeVector direction(1.0,0.0,0.0);
  particleGun->SetParticleMomentumDirection(direction.unit());

  G4double energy = 10.0*MeV;
  particleGun->SetParticleEnergy(energy);

  particleGun->GeneratePrimaryVertex(anEvent);
}
