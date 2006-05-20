#include "StatAccepTestPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"


StatAccepTestPrimaryGeneratorAction::StatAccepTestPrimaryGeneratorAction() {
  G4int n_particle = 1;
  particleGun = new G4ParticleGun( n_particle );

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  particleGun->SetParticleDefinition( particleTable->FindParticle( "geantino" ) );
  particleGun->SetParticleEnergy( 10.0*GeV );

  particleGun->SetParticlePosition( G4ThreeVector( 0.0, 0.0, -2.0*m ) ); //***LOOKHERE***

}


StatAccepTestPrimaryGeneratorAction::~StatAccepTestPrimaryGeneratorAction() {
  delete particleGun;
}


void StatAccepTestPrimaryGeneratorAction::GeneratePrimaries( G4Event* anEvent ) {

  G4ThreeVector v( 0.0, 0.0, 1.0 );  //***LOOKHERE***

  particleGun->SetParticleMomentumDirection( v );
  particleGun->GeneratePrimaryVertex( anEvent );
}


