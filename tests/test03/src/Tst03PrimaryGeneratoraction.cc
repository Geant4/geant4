// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "Tst03PrimaryGeneratorAction.hh"

#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

Tst03PrimaryGeneratorAction::Tst03PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

}

Tst03PrimaryGeneratorAction::~Tst03PrimaryGeneratorAction()
{
  delete particleGun;
}

void Tst03PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
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
