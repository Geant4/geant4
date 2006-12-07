#include "RandomCaloPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "RandomCaloStackingAction.hh"
#include <vector>
#include <string>


// ***LOOKHERE*** These are the particle types that can be used as
//                primary beam particle, on a event-by-event based.
const int RandomCaloPrimaryGeneratorAction::numberCandidateParticles = 17;
const std::string RandomCaloPrimaryGeneratorAction::
nameParticlesVector[ RandomCaloPrimaryGeneratorAction::numberCandidateParticles ] = {
  "mu-", 
  "mu+", 
  "e-",
  "e+", 
  "gamma",
  "pi-", 
  "pi+", 
  "kaon-", 
  "kaon+", 
  "kaon0L", 
  "neutron", 
  "proton",
  "anti_neutron", 
  "anti_proton",
  "deuteron",
  "triton", 
  "alpha"
};


RandomCaloPrimaryGeneratorAction::RandomCaloPrimaryGeneratorAction() :

  infBeamEnergy(   1.0*GeV ),               // ***LOOKHERE***
  supBeamEnergy( 100.0*GeV ),               // ***LOOKHERE***
  // These two constants specify the range of the beam energy.

  rMax( 51.0*cm ),                          // ***LOOKHERE*** 
  zMax( 1346.0*cm ),                        // ***LOOKHERE*** 
  // These two constants specify the cylindrical volume inside
  // which the beam primary particle initial positions can be 
  // randomly drawn.
  // rMax is the radius in x-y plane, whereas zMax is the max
  // value along the z-axis, which is the beam axis.
  // These two constants should be set in between the cylindrical
  // calorimeter and the experimental hall (i.e. the world volume):
  // these values are specified in the class 
  //    RandomCaloDetectorConstruction .

  energyThresholdNoBiasBelow( 10.0*GeV )    // ***LOOKHERE***
  // This constant is the beam energy threshold below which 
  // no biasing is applied for electrons/positrons/gammas/neutrons,
  // in an event-by-event basis.
  // The parameters which specify the biasing, and the application
  // of the biasing itself are in the class RandomCaloStackingAction.

{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun( n_particle );
  particleTable = G4ParticleTable::GetParticleTable();

  G4cout << " --- Info RandomCaloPrimaryGeneratorAction --- " << G4endl
	 << "\t numberCandidateParticles = " 
         << RandomCaloPrimaryGeneratorAction::numberCandidateParticles << G4endl  
	 << "\t candidate primary beam particles: " << G4endl;
  for ( int i = 0; i < RandomCaloPrimaryGeneratorAction::numberCandidateParticles; 
	i++ ) { 
    G4cout << "\t \t" << i << ")  " << nameParticlesVector[ i ] << G4endl;
  }
  G4cout << "\t infBeamEnergy = " << infBeamEnergy / GeV << " GeV" << G4endl
	 << "\t supBeamEnergy = " << supBeamEnergy / GeV << " GeV" << G4endl
	 << "\t rMax = " << rMax / cm << " cm" << G4endl
	 << "\t zMax = " << zMax / cm << " cm" << G4endl
	 << " --------------------------------------------- " << G4endl;

}


RandomCaloPrimaryGeneratorAction::~RandomCaloPrimaryGeneratorAction() {
  delete particleGun;
}


void RandomCaloPrimaryGeneratorAction::GeneratePrimaries( G4Event* anEvent ) {

  // Select randomly the primary particle.
  G4int caseBeamParticle = 
    static_cast< int >( G4UniformRand() * 
			RandomCaloPrimaryGeneratorAction::numberCandidateParticles );
  std::string nameParticle = nameParticlesVector[ caseBeamParticle ];
  particleGun->SetParticleDefinition( particleTable->FindParticle( nameParticle ) );

  // Select randomly the beam energy.
  G4double beamEnergy = infBeamEnergy + 
    G4UniformRand() * ( supBeamEnergy - infBeamEnergy );
  particleGun->SetParticleEnergy( beamEnergy );
  // Do not apply biasing in this event if the beam energy is 
  // below a given threshold.
  if ( beamEnergy < energyThresholdNoBiasBelow ) {
    RandomCaloStackingAction::switchOn( false );    // No biasing in this event.
  } else {  
    RandomCaloStackingAction::switchOn( true );     // Biasing is applied in this event.
  }

  // Select randomly the beam (starting) position
  G4double rPos = G4UniformRand() * rMax;
  G4double phiPos = G4UniformRand() * 2.0 * 3.141592654;
  G4double zPos = G4UniformRand() * zMax;
  if ( G4UniformRand() < 0.5 ) {
    zPos *= -1.0;                // Negative z position.
  }
  G4ThreeVector positionVector( rPos*std::cos( phiPos ), rPos*std::sin( phiPos ), zPos );
  particleGun->SetParticlePosition( positionVector );
  
  // Select randomly the beam direction.
  G4ThreeVector directionVector( G4RandomDirection() );
  particleGun->SetParticleMomentumDirection( directionVector );

  //***DEBUG***
  //G4cout << " RandomCaloPrimaryGeneratorAction::GeneratePrimaries " << G4endl
  //       << "\t particle = " << nameParticle << G4endl
  //	   << "\t beamEnergy = " << beamEnergy / GeV << " GeV" << G4endl
  //	   << "\t positionVector = " << positionVector / cm << " cm" << G4endl
  //       << "\t directionVector = " << directionVector << G4endl;

  particleGun->GeneratePrimaryVertex( anEvent );

}


