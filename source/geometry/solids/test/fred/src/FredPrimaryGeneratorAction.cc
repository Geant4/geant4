//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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
	sprayGun = new SprayParticleGun( );
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
