
#include "DicomPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "DicomGeometry.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

DicomPrimaryGeneratorAction::DicomPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);  	     
}

DicomPrimaryGeneratorAction::~DicomPrimaryGeneratorAction()
{
  delete particleGun;
}

void DicomPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  particleGun->SetParticleDefinition(particle);
  // ---- MGP ---- Numbers in the code should be replaced by const
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.*cm,0.*cm));       
  particleGun->SetParticleEnergy(5.*MeV);
  particleGun->SetParticlePosition(G4ThreeVector(-1.*m,0.,0.));
  particleGun->GeneratePrimaryVertex(anEvent);
}

