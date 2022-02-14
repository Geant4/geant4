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

// (adapted from B1PrimaryGeneratorAction)
// Author: A.Knaian (ara@nklabs.com), N.MacFadden (natemacfadden@gmail.com)

#include "FAPrimaryGeneratorAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
	fParticleGun(0), 
	fWorldBox(0)
{
	G4int n_particle = 1;
	fParticleGun  = new G4ParticleGun(n_particle);

	// default particle kinematic
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle
		= particleTable->FindParticle(particleName="proton");
	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
	fParticleGun->SetParticleEnergy(50.*MeV);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4double worldSizeXY = 0;
	G4double worldSizeZ = 0;

	if (!fWorldBox)
	{
		G4LogicalVolume* worldLV
			= G4LogicalVolumeStore::GetInstance()->GetVolume("World");
		if ( worldLV ) fWorldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
	}

	if ( fWorldBox ) {
		worldSizeXY = fWorldBox->GetXHalfLength()*2.;
		worldSizeZ = fWorldBox->GetZHalfLength()*2.;
	}  
	else  {
		G4ExceptionDescription msg;
		msg << "World volume of box shape not found.\n"; 
		msg << "Perhaps you have changed geometry.\n";
		msg << "The gun will be place at the center.";
		G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
		 "MyCode0002",JustWarning,msg);
	}

	// shoot on XY disk centered on Z-axis behind the cloud
	G4double sigma = worldSizeXY/10.0;			// spread in x and y
	G4double x0 = G4RandGauss::shoot(0,sigma);
	G4double y0 = G4RandGauss::shoot(0,sigma);
	G4double z0 = 0.95 * (-0.5) * worldSizeZ;
	
	fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

	fParticleGun->GeneratePrimaryVertex(anEvent);
}

