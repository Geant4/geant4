//
// FredPrimaryGeneratorAction.cc
//
// Implementation of Fred's generator
//

#include "FredPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "SprayParticleGun.hh"
#include "GridParticleGun.hh"
#include "G4ParticleGun.hh"
#include "G4UImanager.hh"
#include "globals.hh"

//
// (Constructor)
//
FredPrimaryGeneratorAction::FredPrimaryGeneratorAction( FredMessenger *ourMessenger )
{
	//
	// Fred has a few interesting ways of generating particles
	//
	sprayGun = new SprayParticleGun( 1 );
	gridGun  = new GridParticleGun( );
	g4Gun    = new G4ParticleGun();
	
	//
	// We ask fred's main messenger which to use
	//
	messenger = ourMessenger;
}

//
// (Destructor)
//
FredPrimaryGeneratorAction::~FredPrimaryGeneratorAction()
{
	delete sprayGun;
	delete gridGun;
	delete g4Gun;
}

//
// GeneratePrimaries
//
void FredPrimaryGeneratorAction::GeneratePrimaries( G4Event *anEvent )
{
	switch (messenger->SelectedGun()) {
		case SPRAY:
		sprayGun->GeneratePrimaryVertex( anEvent );
		break;
		
		case GRID:
		gridGun->GeneratePrimaryVertex( anEvent );
		break;
		
		case G4:
		g4Gun->GeneratePrimaryVertex( anEvent );
		break;
	}
}
