#include "MyPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"


MyPrimaryGeneratorAction::MyPrimaryGeneratorAction() {
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  particleGun->SetParticleDefinition( particleTable->FindParticle("geantino") );
  particleGun->SetParticleEnergy( 10.0*GeV );
  particleGun->SetParticlePosition( G4ThreeVector(-1.5*m, 0.0, 0.0) );
}


MyPrimaryGeneratorAction::~MyPrimaryGeneratorAction() {
  delete particleGun;
}


void MyPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  G4ThreeVector v(1.0, 0.0, 0.0);
  particleGun->SetParticleMomentumDirection(v);
  particleGun->GeneratePrimaryVertex(anEvent);
}


