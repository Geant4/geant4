#include "StatAccepTestPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"


StatAccepTestPrimaryGeneratorAction::StatAccepTestPrimaryGeneratorAction() {
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  particleGun->SetParticleDefinition( particleTable->FindParticle("geantino") );
  particleGun->SetParticleEnergy( 10.0*GeV );
  particleGun->SetParticlePosition( G4ThreeVector(-2.0*m, 0.0, 0.0) );
}


StatAccepTestPrimaryGeneratorAction::~StatAccepTestPrimaryGeneratorAction() {
  delete particleGun;
}


void StatAccepTestPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  G4ThreeVector v(1.0, 0.0, 0.0);
  particleGun->SetParticleMomentumDirection(v);
  particleGun->GeneratePrimaryVertex(anEvent);
}


