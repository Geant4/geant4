#include "MyPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include <vector>
#include <string>


// ***LOOKHERE*** These are the particle types that can be used as
//                primary beam particle, on a event-by-event based.
const int MyPrimaryGeneratorAction::numberCandidateParticles = 29;
const std::string MyPrimaryGeneratorAction::
nameParticlesVector[ MyPrimaryGeneratorAction::numberCandidateParticles ] = {
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
  "alpha",
  "lambda",
  "sigma+",
  "sigma-",
  "xi-",
  "xi0",
  "anti_lambda",
  "anti_sigma+",
  "anti_sigma-",
  "anti_xi-",
  "anti_xi0",
  "omega-",
  "anti_omega-"
};


MyPrimaryGeneratorAction::MyPrimaryGeneratorAction() :
  // These two constants specify the range of the beam energy.
  infBeamEnergy(   1.0*GeV ),               // ***LOOKHERE***
  supBeamEnergy( 100.0*GeV )                // ***LOOKHERE***
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun( n_particle );
  particleTable = G4ParticleTable::GetParticleTable();

  G4cout << " --- Info MyPrimaryGeneratorAction --- " << G4endl
	 << "\t numberCandidateParticles = " 
         << MyPrimaryGeneratorAction::numberCandidateParticles << G4endl  
	 << "\t candidate primary beam particles: " << G4endl;
  for ( int i = 0; i < MyPrimaryGeneratorAction::numberCandidateParticles; 
	i++ ) { 
    G4cout << "\t \t" << i << ")  " << nameParticlesVector[ i ] << G4endl;
  }
  G4cout << "\t infBeamEnergy = " << infBeamEnergy / GeV << " GeV" << G4endl
	 << "\t supBeamEnergy = " << supBeamEnergy / GeV << " GeV" << G4endl
	 << " --------------------------------------------- " << G4endl;
}


MyPrimaryGeneratorAction::~MyPrimaryGeneratorAction() {
  delete particleGun;
}


void MyPrimaryGeneratorAction::GeneratePrimaries( G4Event* anEvent ) {

  // Select randomly the primary particle.
  G4int caseBeamParticle = 
    static_cast< int >( G4UniformRand() * 
			MyPrimaryGeneratorAction::numberCandidateParticles );
  std::string nameParticle = nameParticlesVector[ caseBeamParticle ];

  // Select randomly the beam energy.
  G4double beamEnergy = infBeamEnergy + 
    G4UniformRand() * ( supBeamEnergy - infBeamEnergy );

  // Beam position: always the origin.
  G4ThreeVector positionVector( 0.0, 0.0, 0.0 );
  
  // Select randomly the beam direction.
  G4ThreeVector directionVector( G4RandomDirection() );


  //***LOOKHERE*** If you want to specify a given beam particle,
  //               with a given beam energy, position, and direction:
  //
  //nameParticle = "geantino"
  //nameParticle = "chargedgeantino"
  //nameParticle = "mu-"
  //nameParticle = "mu+"
  //nameParticle = "e-"
  //nameParticle = "e+"
  //nameParticle = "gamma"
  //nameParticle = "pi-";
  //nameParticle = "pi+"
  //nameParticle = "kaon-"
  //nameParticle = "kaon+"
  //nameParticle = "kaon0L"
  //nameParticle = "neutron"
  //nameParticle = "proton"
  //
  //beamEnergy = 20.0*GeV;
  //  
  //directionVector = G4ThreeVector( 0.0, 1.0, 0.0 );
  //***endLOOKHERE***


  particleGun->SetParticleDefinition( particleTable->FindParticle( nameParticle ) );
  particleGun->SetParticleEnergy( beamEnergy );
  particleGun->SetParticlePosition( positionVector );
  particleGun->SetParticleMomentumDirection( directionVector );

  //***DEBUG***
  G4cout << " MyPrimaryGeneratorAction::GeneratePrimaries " << G4endl
         << "\t particle = " << nameParticle << G4endl
         << "\t beamEnergy = " << beamEnergy / GeV << " GeV" << G4endl
	 << "\t positionVector = " << positionVector / cm << " cm" << G4endl
         << "\t directionVector = " << directionVector << G4endl;

  particleGun->GeneratePrimaryVertex( anEvent );

}


