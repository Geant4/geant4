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

#include "G4INCLNDeltaToNSKChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NDeltaToNSKChannel::angularSlope = 2.;
	
	NDeltaToNSKChannel::NDeltaToNSKChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NDeltaToNSKChannel::~NDeltaToNSKChannel(){}
	
	void NDeltaToNSKChannel::fillFinalState(FinalState *fs) {
        // D++ p -> p S+ K+ (6)
        //
        // D++ n -> p S+ K0 (3)
        // D++ n -> p S0 K+ (3)
        // D++ n -> n S+ K+ (3)
        //
        // D+  p -> p S+ K0 (2)
        // D+  p -> p S0 K+ (2)
        // D+  p -> n S+ K+ (3)
        //
        // D+  n -> p S0 K0 (3)
        // D+  n -> p S- K+ (2)
        // D+  n -> n S+ K0 (2)
        // D+  n -> n S0 K+ (2)
        
        
        Particle *delta;
        
        if (particle1->isResonance()) {
            delta = particle1;
        }
        else {
            delta = particle2;
        }
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(particle1, particle2);
		
		const G4int iso = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
		const G4int iso_d = ParticleTable::getIsospin(delta->getType());
		
		ParticleType KaonType;
		ParticleType NucleonType;
		ParticleType SigmaType;
		
		const G4double rdm = Random::shoot();
		
		if(std::abs(iso) == 4){// D++ p
			KaonType = ParticleTable::getKaonType(iso/4);
			NucleonType = ParticleTable::getNucleonType(iso/4);
			SigmaType = ParticleTable::getSigmaType(iso/2);
		}
		else if(iso == 0){// D+  n
			if(rdm*9 < 3){
				KaonType = ParticleTable::getKaonType(-iso_d);
				NucleonType = ParticleTable::getNucleonType(iso_d);
				SigmaType = SigmaZero;
			}
			else if(rdm*9 < 5){
				KaonType = ParticleTable::getKaonType(iso_d);
				NucleonType = ParticleTable::getNucleonType(iso_d);
				SigmaType = ParticleTable::getSigmaType(-2*iso_d);
			}
			else if(rdm*9 < 7){
				KaonType = ParticleTable::getKaonType(-iso_d);
				NucleonType = ParticleTable::getNucleonType(-iso_d);
				SigmaType = ParticleTable::getSigmaType(2*iso_d);
			}
			else{
				KaonType = ParticleTable::getKaonType(iso_d);
				NucleonType = ParticleTable::getNucleonType(-iso_d);
				SigmaType = SigmaZero;
			}
		}
		else if(ParticleTable::getIsospin(particle1->getType()) == ParticleTable::getIsospin(particle2->getType())){// D+  p
			if(rdm*7 < 2){
				KaonType = ParticleTable::getKaonType(-iso/2);
				NucleonType = ParticleTable::getNucleonType(iso/2);
				SigmaType = ParticleTable::getSigmaType(iso);
			}
			else if(rdm*7 < 4){
				KaonType = ParticleTable::getKaonType(iso/2);
				NucleonType = ParticleTable::getNucleonType(iso/2);
				SigmaType = SigmaZero;
			}
			else{
				KaonType = ParticleTable::getKaonType(iso/2);
				NucleonType = ParticleTable::getNucleonType(-iso/2);
				SigmaType = ParticleTable::getSigmaType(iso);
			}
		}
		else{// D++ n 
			if(rdm*3 < 1){
				KaonType = ParticleTable::getKaonType(-iso/2);
				NucleonType = ParticleTable::getNucleonType(iso/2);
				SigmaType = ParticleTable::getSigmaType(iso);
			}
			else if(rdm*3 < 2){
				KaonType = ParticleTable::getKaonType(iso/2);
				NucleonType = ParticleTable::getNucleonType(iso/2);
				SigmaType = SigmaZero;
			}
			else{
				KaonType = ParticleTable::getKaonType(iso/2);
				NucleonType = ParticleTable::getNucleonType(-iso/2);
				SigmaType = ParticleTable::getSigmaType(iso);
			}
		}
		
		particle1->setType(NucleonType);
		particle2->setType(SigmaType);
		
		ParticleList list;
		list.push_back(particle1);
		list.push_back(particle2);
		const ThreeVector &rcol = particle2->getPosition();
		const ThreeVector zero;
		Particle *kaon = new Particle(KaonType,zero,rcol);
		list.push_back(kaon);
		
		PhaseSpaceGenerator::generateBiased(sqrtS, list, 0, angularSlope);
		
		fs->addModifiedParticle(particle1);
		fs->addModifiedParticle(particle2);
		fs->addCreatedParticle(kaon);
		
	}
}
