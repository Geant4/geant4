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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLNDeltaToNNKKbChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NDeltaToNNKKbChannel::angularSlope = 2.;
	
	NDeltaToNNKKbChannel::NDeltaToNNKKbChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NDeltaToNNKKbChannel::~NDeltaToNNKKbChannel(){}
	
	void NDeltaToNNKKbChannel::fillFinalState(FinalState *fs) {
		
        // ratio
        // D++ p -> p p K+ K0b	(6)
        //
        // D++ n -> p p K+ K-	(3)
        // D++ n -> p p K0 K0b	(3)
        // D++ n -> p n K+ K0b	(3)
        //
        // D+  p -> p p K+ K-	(3)
        // D+  p -> p p K0 K0b	(1)
        // D+  p -> p n K+ K0b	(3)
        //
        // D+  n -> p p K0 K-	(2)
        // D+  n -> p n K+ K-	(1)
        // D+  n -> p n K0 K0b	(3)
        // D+  n -> n n K+ K0b	(2)
        //
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(particle1, particle2);
		
		const G4int iso = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
		G4int iso_n;
		if(particle1->isNucleon()) iso_n = ParticleTable::getIsospin(particle1->getType());
		else iso_n = ParticleTable::getIsospin(particle2->getType());
		
		ParticleType Nucleon1Type;
		ParticleType Nucleon2Type;
		ParticleType KaonType;
		ParticleType antiKaonType;
		
		const G4double rdm = Random::shoot();
		
		if(std::abs(iso) == 4){// D++ p
			Nucleon1Type = ParticleTable::getNucleonType(iso/4);
			Nucleon2Type = ParticleTable::getNucleonType(iso/4);
			KaonType = ParticleTable::getKaonType(iso/4);
			antiKaonType = ParticleTable::getAntiKaonType(iso/4);
		}
		else if(iso == 0){// D+  n
			if(rdm*8 < 2){
				Nucleon1Type = Proton;
				Nucleon2Type = Proton;
				KaonType = KZero;
				antiKaonType = KMinus;
			}
			else if(rdm*8 < 3){
				Nucleon1Type = Proton;
				Nucleon2Type = Neutron;
				KaonType = ParticleTable::getKaonType(-iso_n);
				antiKaonType = ParticleTable::getAntiKaonType(iso_n);
			}
			else if(rdm*8 < 6){
				Nucleon1Type = Proton;
				Nucleon2Type = Neutron;
				KaonType = ParticleTable::getKaonType(iso_n);
				antiKaonType = ParticleTable::getAntiKaonType(-iso_n);
			}
			else{
				Nucleon1Type = Neutron;
				Nucleon2Type = Neutron;
				KaonType = KPlus;
				antiKaonType = KZeroBar;
			}
		}
		else if(ParticleTable::getIsospin(particle1->getType()) == ParticleTable::getIsospin(particle2->getType())){// D+  p
			if(rdm*3 < 1){
				Nucleon1Type = ParticleTable::getNucleonType(iso/2);
				Nucleon2Type = ParticleTable::getNucleonType(iso/2);
				KaonType = ParticleTable::getKaonType(iso/2);
				antiKaonType = ParticleTable::getAntiKaonType(-iso/2);
			}
			else if(rdm*3 < 2){
				Nucleon1Type = ParticleTable::getNucleonType(iso/2);
				Nucleon2Type = ParticleTable::getNucleonType(iso/2);
				KaonType = ParticleTable::getKaonType(-iso/2);
				antiKaonType = ParticleTable::getAntiKaonType(iso/2);
			}
			else{
				Nucleon1Type = ParticleTable::getNucleonType(iso/2);
				Nucleon2Type = ParticleTable::getNucleonType(-iso/2);
				KaonType = ParticleTable::getKaonType(iso/2);
				antiKaonType = ParticleTable::getAntiKaonType(iso/2);
			}
		}
		else{// D++ n 
			if(rdm*5 < 2){
				Nucleon1Type = ParticleTable::getNucleonType(iso/2);
				Nucleon2Type = ParticleTable::getNucleonType(iso/2);
				KaonType = ParticleTable::getKaonType(iso/2);
				antiKaonType = ParticleTable::getAntiKaonType(-iso/2);
			}
			else if(rdm*5 < 4){
				Nucleon1Type = ParticleTable::getNucleonType(iso/2);
				Nucleon2Type = ParticleTable::getNucleonType(iso/2);
				KaonType = ParticleTable::getKaonType(-iso/2);
				antiKaonType = ParticleTable::getAntiKaonType(iso/2);
			}
			else{
				Nucleon1Type = ParticleTable::getNucleonType(iso/2);
				Nucleon2Type = ParticleTable::getNucleonType(-iso/2);
				KaonType = ParticleTable::getKaonType(iso/2);
				antiKaonType = ParticleTable::getAntiKaonType(iso/2);
			}
		}
		
				
		particle1->setType(Nucleon1Type);
		particle2->setType(Nucleon2Type);
		
		ParticleList list;
		list.push_back(particle1);
		list.push_back(particle2);
		const ThreeVector &rcol1 = particle1->getPosition();
		const ThreeVector zero1;
		const ThreeVector &rcol2 = particle2->getPosition();
		const ThreeVector zero2;
		Particle *kaon = new Particle(KaonType,zero1,rcol1);
		Particle *antikaon = new Particle(antiKaonType,zero2,rcol2);
		list.push_back(kaon);
		list.push_back(antikaon);
		
		PhaseSpaceGenerator::generateBiased(sqrtS, list, 0, angularSlope);
		
		fs->addModifiedParticle(particle1);
		fs->addModifiedParticle(particle2);
		fs->addCreatedParticle(kaon);
		fs->addCreatedParticle(antikaon);
	
	}
}
