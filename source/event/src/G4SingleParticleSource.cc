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

///////////////////////////////////////////////////////////////////////////////
//
// MODULE:        G4SingleParticleSource.hh
//
// Version:      1.0
// Date:         5/02/04
// Author:       Fan Lei 
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// Version 1.0, 05/02/2004, Fan Lei, Created.
//      Based on the G4GeneralParticleSource class in Geant4 v6.0
//
///////////////////////////////////////////////////////////////////////////////
//
#include <cmath>

#include "G4SingleParticleSource.hh"

#include "G4SystemOfUnits.hh"
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"
#include "G4Geantino.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4AutoLock.hh"

G4SingleParticleSource::part_prop_t::part_prop_t() {
  //definition = G4Geantino::GeantinoDefinition();
  momentum_direction = G4ParticleMomentum(1,0,0);
  energy = 1.*MeV;
  position = G4ThreeVector();
}

G4SingleParticleSource::G4SingleParticleSource() {
//	// Initialise all variables
//	// Position distribution Variables
//
	NumberOfParticlesToBeGenerated = 1;
	definition = G4Geantino::GeantinoDefinition();
//	G4ThreeVector zero;
//	particle_momentum_direction = G4ParticleMomentum(1, 0, 0);
//	particle_energy = 1.0 * MeV;
//	particle_position = zero;
//	particle_time = 0.0;
//	particle_polarization = zero;
	charge = 0.0;
	time = 0;
	polarization = G4ThreeVector();

	biasRndm = new G4SPSRandomGenerator();
	posGenerator = new G4SPSPosDistribution();
	posGenerator->SetBiasRndm(biasRndm);
	angGenerator = new G4SPSAngDistribution();
	angGenerator->SetPosDistribution(posGenerator);
	angGenerator->SetBiasRndm(biasRndm);
	eneGenerator = new G4SPSEneDistribution();
	eneGenerator->SetBiasRndm(biasRndm);

	// verbosity
	verbosityLevel = 0;
    
    G4MUTEXINIT(mutex);

}

G4SingleParticleSource::~G4SingleParticleSource() {
	delete biasRndm;
	delete posGenerator;
	delete angGenerator;
	delete eneGenerator;
    G4MUTEXDESTROY(mutex);
}

void G4SingleParticleSource::SetVerbosity(int vL) {
        G4AutoLock l(&mutex);
	verbosityLevel = vL;
	posGenerator->SetVerbosity(vL);
	angGenerator->SetVerbosity(vL);
	eneGenerator->SetVerbosity(vL);
	//G4cout << "Verbosity Set to: " << verbosityLevel << G4endl;
}

void G4SingleParticleSource::SetParticleDefinition(
		G4ParticleDefinition* aParticleDefinition) {
	definition = aParticleDefinition;
	charge = aParticleDefinition->GetPDGCharge();
}

void G4SingleParticleSource::GeneratePrimaryVertex(G4Event *evt) {

    //G4AutoLock l(&mutex);
        //part_prop_t& pp = ParticleProperties.Get();
	if (definition == NULL) {
		//TODO: Should this rise an exception???
		return ;
	}
		//return;

	if (verbosityLevel > 1)
		G4cout << " NumberOfParticlesToBeGenerated: "
				<<NumberOfParticlesToBeGenerated << G4endl;

        part_prop_t& pp = ParticleProperties.Get();
	// Position stuff
	pp.position = posGenerator->GenerateOne();

	// create a new vertex
	G4PrimaryVertex* vertex = new G4PrimaryVertex(pp.position,time);

	for (G4int i = 0; i < NumberOfParticlesToBeGenerated; i++) {
		// Angular stuff
		pp.momentum_direction = angGenerator->GenerateOne();
		// Energy stuff
		pp.energy = eneGenerator->GenerateOne(definition);

		if (verbosityLevel >= 2)
			G4cout << "Creating primaries and assigning to vertex" << G4endl;
		// create new primaries and set them to the vertex
		G4double mass = definition->GetPDGMass();
		G4PrimaryParticle* particle =
		  new G4PrimaryParticle(definition);
		particle->SetKineticEnergy(pp.energy );
		particle->SetMass( mass );
		particle->SetMomentumDirection( pp.momentum_direction );
		particle->SetCharge( charge );
		particle->SetPolarization(polarization.x(),
					  polarization.y(),
					  polarization.z());
		if (verbosityLevel > 1) {
			G4cout << "Particle name: "
					<< definition->GetParticleName() << G4endl;
			G4cout << "       Energy: " << pp.energy << G4endl;
			G4cout << "     Position: " << pp.position << G4endl;
			G4cout << "    Direction: " << pp.momentum_direction
					<< G4endl;
		}
		// Set bweight equal to the multiple of all non-zero weights
		G4double weight = eneGenerator->GetWeight()*biasRndm->GetBiasWeight();
		// pass it to primary particle
		particle->SetWeight(weight);

		vertex->SetPrimary(particle);

	}
	// now pass the weight to the primary vertex. CANNOT be used here!
	//  vertex->SetWeight(particle_weight);
	evt->AddPrimaryVertex(vertex);
	if (verbosityLevel > 1)
		G4cout << " Primary Vetex generated !" << G4endl;
}

