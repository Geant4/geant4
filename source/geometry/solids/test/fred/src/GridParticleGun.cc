//
// GridParticleGun.cc
//
// Implementation of grid particle gun
//

#include "GridParticleGun.hh"
#include "GridParticleGunMessenger.hh"

#include "G4PrimaryParticle.hh"
#include "G4ParticleTypes.hh"
#include "G4Event.hh"

//
// Constructor (no arguments)
//
GridParticleGun::GridParticleGun()
{
	SetDefaults();
}

//
// Constructor
//
GridParticleGun::GridParticleGun( G4ThreeVector theDirection,
				  G4ThreeVector theOrigin,
				  G4ThreeVector theGrid1,
				  G4ThreeVector theGrid2,
				  G4int theN1, G4int theN2 )
{
	SetDefaults();
	
	direction = theDirection;
	origin    = theOrigin;
	grid1	  = theGrid1;
	grid2     = theGrid2;
	n1        = theN1;
	n2        = theN2;
	
	if (n1*n2 == 0) 
		G4Exception( "GridParticleGun created with n1 or n2 zero" );
}
	
	
//
// Destructor
//
GridParticleGun::~GridParticleGun() {
	delete messenger;
}


//
// SetDefaults
//
// Set defaults and create associated messenger
//
void GridParticleGun::SetDefaults()
{
	direction = G4ThreeVector(0,0,1);
	origin = G4ThreeVector(-2*m,-2*m,-2*m);
	grid1 = G4ThreeVector(4*m,0,0);
	grid2 = G4ThreeVector(0,4*m,0);
	n1 = 10;
	n2 = 10;
	
	
	particleType = G4Geantino::GeantinoDefinition();
	
	messenger = new GridParticleGunMessenger( this );
}


//
// GeneratePrimaryVertex
//
// Make the event.
//
// Note the strange notation of G4Event: we are allowed to have as many
// "primary" vertices as we wish.
//
void GridParticleGun::GeneratePrimaryVertex( G4Event *evt )
{
	G4ThreeVector 	delta1 = grid1*(1.0/((G4double) n1)),
					delta2 = grid2*(1.0/((G4double) n2));

	G4int i = 0;
	G4ThreeVector iStart = origin;
	do {
		G4int j = 0;
		G4ThreeVector jStart = iStart;
		do {
		
			//
			// New particle: make it's vertex
			//
			G4PrimaryVertex *vertex = new G4PrimaryVertex( jStart, 0.0 );
			
			//
			// ...then make it's self
			//
			G4PrimaryParticle *particle = new G4PrimaryParticle( particleType, direction.x(), direction.y(), direction.z() );
			
			//
			// ...add it to the vertex, and add the vertex to the event
			//
			vertex->SetPrimary( particle );
			evt->AddPrimaryVertex( vertex );
			
		} while( jStart += delta2, j++ < n2 );
	} while( iStart += delta1, i++ < n1 );
}
			
