//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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
