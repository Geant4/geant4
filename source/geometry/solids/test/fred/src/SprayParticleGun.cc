//
// SprayParticleGun
//
// Implementation of test particle gun
//

#include "SprayParticleGun.hh"
#include "SprayParticleGunMessenger.hh"

#include "G4PrimaryParticle.hh"
#include "G4ParticleTypes.hh"
#include "G4Event.hh"

//
// Constructor (no arguments)
//
SprayParticleGun::SprayParticleGun()
{
	SetDefaults();
}

//
// Constructor (position specified)
//
SprayParticleGun::SprayParticleGun( G4ThreeVector inStart )
{
	SetDefaults();
	start = inStart;
}

//
// Destructor
//
SprayParticleGun::~SprayParticleGun() {
	delete messenger;
}

//
// SetDefaults
//
// Set default parameters
//
void SprayParticleGun::SetDefaults()
{
	start = G4ThreeVector(0,0,0);
	spray[0] = 7;
	spray[1] = 7;
	spray[2] = 2;
	
	particleType = G4Geantino::GeantinoDefinition();
	
  	messenger = new SprayParticleGunMessenger(this);
}

//
// GeneratePrimaryVertex
//
// Make the event
//
void SprayParticleGun::GeneratePrimaryVertex( G4Event *evt )
{
	//
	// Create a new vertex
	//
	G4PrimaryVertex *vertex = new G4PrimaryVertex( start, 0.0 );
	
	//
	// Loop over the three dimensions, creating particles
	// as specified
	//
	G4double px = -1.0;
	G4int    ix = 1;
	do {
		if (!(spray[0]&ix)) continue;
		
		G4double py = -1.0;
		G4int    iy = 1;
		do {
			if (!(spray[1]&iy)) continue;
		
			G4double pz = -1.0;
			G4int    iz = 1; 
			do {
				if (!(spray[2]&iz)) continue;
			
				if (ix==2 && iy==2 && iz==2) continue;
				
				G4PrimaryParticle *particle = new G4PrimaryParticle( particleType, px, py, pz );
				vertex->SetPrimary( particle );
			} while( pz += 1.0, iz <<= 1, iz < 8 );
		} while( py += 1.0, iy <<= 1, iy < 8 );
	} while( px += 1.0, ix <<= 1, ix < 8 );

	//
	// Give this sucker to the event
	//
	evt->AddPrimaryVertex( vertex );
}
