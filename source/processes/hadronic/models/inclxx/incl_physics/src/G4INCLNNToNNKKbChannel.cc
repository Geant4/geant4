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

#include "G4INCLNNToNNKKbChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NNToNNKKbChannel::angularSlope = 2.; // What is the exact effect? Sould be check
	
	NNToNNKKbChannel::NNToNNKKbChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NNToNNKKbChannel::~NNToNNKKbChannel(){}
	
	void NNToNNKKbChannel::fillFinalState(FinalState *fs) {
		
        // pp -> pp K+ K- (1)
        // pp -> pp K0 K0 (1)
        // pp -> pn K+ K0 (4)
        // pn -> pp K0 K- (4)
        // pn -> pn K+ K- (9)
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(particle1, particle2);
		
		const G4int iso = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
		
		ParticleType Kaon1Type;
		ParticleType Kaon2Type;
		G4double rdm = Random::shoot();
		
		if(iso == 2){
			if(rdm*6. < 1.){
				Kaon1Type = KPlus;
				Kaon2Type = KMinus;
			}
			else if(rdm*6. < 2){
				Kaon1Type = KZero;
				Kaon2Type = KZeroBar;
			}
			else{
				Kaon1Type = KPlus;
				Kaon2Type = KZeroBar;
				particle1->setType(Neutron);
			}
		}
		else if(iso == -2){
			if(rdm*6. < 1.){
				Kaon1Type = KPlus;
				Kaon2Type = KMinus;
			}
			else if(rdm*6. < 2){
				Kaon1Type = KZero;
				Kaon2Type = KZeroBar;
			}
			else{
				Kaon1Type = KZero;
				Kaon2Type = KMinus;
				particle1->setType(Proton);
			}
		}
		else if(rdm*26. < 9.){
			Kaon1Type = KPlus;
			Kaon2Type = KMinus;
		}
		else if(rdm*26. < 18.){
			Kaon1Type = KZero;
			Kaon2Type = KZeroBar;
		}
		else if(rdm*26. < 22.){
			Kaon1Type = KZero;
			Kaon2Type = KMinus;
			particle1->setType(Proton);
			particle2->setType(Proton);
		}
		else{
			Kaon1Type = KPlus;
			Kaon2Type = KZeroBar;
			particle1->setType(Neutron);
			particle2->setType(Neutron);
		}
		
		ParticleList list;
		list.push_back(particle1);
		list.push_back(particle2);
		const ThreeVector &rcol1 = particle1->getPosition();
		const ThreeVector &rcol2 = particle2->getPosition();
		const ThreeVector zero;
		Particle *kaon = new Particle(Kaon1Type,zero,rcol1);
		Particle *antikaon = new Particle(Kaon2Type,zero,rcol2);
		list.push_back(kaon);
		list.push_back(antikaon);
		
		PhaseSpaceGenerator::generateBiased(sqrtS, list, 0, angularSlope);
		
		INCL_DEBUG("NNToNNKKb " << (kaon->getMomentum().theta()) * 180. / G4INCL::Math::pi << '\n');
		
		fs->addModifiedParticle(particle1);
		fs->addModifiedParticle(particle2);
		fs->addCreatedParticle(kaon);
		fs->addCreatedParticle(antikaon);
		
	}
}
