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

#include "G4INCLNDeltaToNLKChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NDeltaToNLKChannel::angularSlope = 2.;
	
	NDeltaToNLKChannel::NDeltaToNLKChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NDeltaToNLKChannel::~NDeltaToNLKChannel(){}
	
	void NDeltaToNLKChannel::fillFinalState(FinalState *fs) {
        // D++ n -> p L K+ (3)
        //
        // D+  p -> p L K+ (1)
        //
        // D+  n -> p L K0 (1)
        // D+  n -> n L K+ (1)
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(particle1, particle2);
		
		const G4int iso = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
		
		ParticleType KaonType;
		ParticleType NucleonType;
		particle2->setType(Lambda);
		
		if(std::abs(iso) == 2){// D++ n and D+ p
			KaonType = ParticleTable::getKaonType(iso/2);
			NucleonType = ParticleTable::getNucleonType(iso/2);
		}
		else if(Random::shoot() < 0.5){// D+  n
			KaonType = KPlus;
			NucleonType = Neutron;
		}
		else{// D+  n
			KaonType = KZero;
			NucleonType = Proton;
		}
		particle1->setType(NucleonType);
		
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
