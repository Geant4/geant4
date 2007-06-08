#include "DetectorSlicePrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"


DetectorSlicePrimaryGeneratorAction::DetectorSlicePrimaryGeneratorAction() {
  G4int n_particle = 1;
  particleGun = new G4ParticleGun( n_particle );

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  particleGun->SetParticleDefinition( particleTable->FindParticle( "geantino" ) );
  particleGun->SetParticleEnergy( 10.0*GeV );

  particleGun->SetParticlePosition( G4ThreeVector( 0.0, 0.0, -1.0*mm ) ); //***LOOKHERE***

}


DetectorSlicePrimaryGeneratorAction::~DetectorSlicePrimaryGeneratorAction() {
  delete particleGun;
}


void DetectorSlicePrimaryGeneratorAction::GeneratePrimaries( G4Event* anEvent ) {

  G4ThreeVector v( 0.0, 0.0, 1.0 );  //***LOOKHERE***

  particleGun->SetParticleMomentumDirection( v );
  particleGun->GeneratePrimaryVertex( anEvent );
}


