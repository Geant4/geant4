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

#include "G4INCLNDeltaToDeltaLKChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NDeltaToDeltaLKChannel::angularSlope = 2.;
	
	NDeltaToDeltaLKChannel::NDeltaToDeltaLKChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NDeltaToDeltaLKChannel::~NDeltaToDeltaLKChannel(){}
	
	G4double NDeltaToDeltaLKChannel::sampleDeltaMass(G4double ecm) {
		const G4double maxDeltaMass = ecm - ParticleTable::effectiveLambdaMass - ParticleTable::effectiveKaonMass - 1.0;
		const G4double maxDeltaMassRndm = std::atan((maxDeltaMass-ParticleTable::effectiveDeltaMass)*2./ParticleTable::effectiveDeltaWidth);
		const G4double deltaMassRndmRange = maxDeltaMassRndm - ParticleTable::minDeltaMassRndm;
// assert(deltaMassRndmRange>0.);

		G4double y=ecm*ecm;
		G4double q2=(y-1.157776E6)*(y-6.4E5)/y/4.0; // 1.157776E6 = 1076^2 = (mNucleon + mPion)^2, 6.4E5 = 800^2 = (mNucleon - mPion)^2
		G4double q3=std::pow(std::sqrt(q2), 3.);
		const G4double f3max=q3/(q3+5.832E6); // 5.832E6 = 180^3 = ???^3
		G4double x;

		G4int nTries = 0;
		G4bool success = false;
		while(!success) { /* Loop checking, 10.07.2015, D.Mancusi */
			if(++nTries >= 100000) {
				INCL_WARN("NDeltaToDeltaLKChannel::sampleDeltaMass loop was stopped because maximum number of tries was reached. Minimum delta mass "
						  << ParticleTable::minDeltaMass << " MeV with CM energy " << ecm << " MeV may be unphysical." << '\n');
				return ParticleTable::minDeltaMass;
			}
			
			G4double rndm = ParticleTable::minDeltaMassRndm + Random::shoot() * deltaMassRndmRange;
			y = std::tan(rndm);
			x = ParticleTable::effectiveDeltaMass + 0.5*ParticleTable::effectiveDeltaWidth*y;
// assert(x>=ParticleTable::minDeltaMass && ecm >= x + ParticleTable::effectiveLambdaMass + ParticleTable::effectiveKaonMass + 1.0);
			
			// generation of the delta mass with the penetration factor
			// (see prc56(1997)2431)
			y=x*x;
			q2=(y-1.157776E6)*(y-6.4E5)/y/4.0; // 1.157776E6 = 1076^2 = (mNucleon + mPion)^2, 6.4E5 = 800^2 = (mNucleon - mPion)^2
			q3=std::pow(std::sqrt(q2), 3.);
			const G4double f3=q3/(q3+5.832E6); // 5.832E6 = 180^3 = ???^3
			rndm = Random::shoot();
			if (rndm*f3max < f3)
			success = true;
		}
		return x;
	}
	
	void NDeltaToDeltaLKChannel::fillFinalState(FinalState *fs) {
        // D++ p -> L K+ D++ (4)
        //
        // D++ n -> L K+ D+  (3)
        // D++ n -> L K0 D++ (4)
        //
        // D+  p -> L K0 D++ (3)
        // D+  p -> L K+ D+  (2)
        //
        // D+  n -> L K+ D0  (4)
        // D+  n -> L K0 D+  (2)
        
        Particle *delta;
        Particle *nucleon;
        
        if (particle1->isResonance()) {
            delta = particle1;
            nucleon = particle2;
        }
        else {
            delta = particle2;
            nucleon = particle1;
        }
		
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(particle1, particle2);
		
		const G4int iso = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
		const G4int iso_d = ParticleTable::getIsospin(delta->getType());
		const G4double rdm = Random::shoot();
		
/*		const G4double m1 = particle1->getMass();
		const G4double m2 = particle2->getMass();
		const G4double pLab = KinematicsUtils::momentumInLab(particle1, particle2);*/
		
		ParticleType KaonType;
		ParticleType DeltaType;
		nucleon->setType(Lambda);
		
		if(std::abs(iso) == 4){// D++ p
			KaonType = ParticleTable::getKaonType(iso/4);
			DeltaType = ParticleTable::getDeltaType(3*iso/4);
		}
		else if(iso == 0){// D+  n
			if(rdm*3 < 2){
				KaonType = ParticleTable::getKaonType(iso_d);
				DeltaType = ParticleTable::getDeltaType(-iso_d);
			}
			else{
				KaonType = ParticleTable::getKaonType(-iso_d);
				DeltaType = ParticleTable::getDeltaType(iso_d);
			}
		}
		else if(ParticleTable::getIsospin(particle1->getType()) == ParticleTable::getIsospin(particle2->getType())){// D+  p
			if(rdm*5 < 3){
				KaonType = ParticleTable::getKaonType(-iso/2);
				DeltaType = ParticleTable::getDeltaType(3*iso/2);
			}
			else{
				KaonType = ParticleTable::getKaonType(iso/2);
				DeltaType = ParticleTable::getDeltaType(iso/2);
			}
		}
		else{// D++ n 
			if(rdm*7 < 3){
				KaonType = ParticleTable::getKaonType(iso/2);
				DeltaType = ParticleTable::getDeltaType(iso/2);
			}
			else{
				KaonType = ParticleTable::getKaonType(-iso/2);
				DeltaType = ParticleTable::getDeltaType(3*iso/2);
			}
		}
		
		delta->setType(DeltaType);
		delta->setMass(sampleDeltaMass(sqrtS));
		
		ParticleList list;
		list.push_back(delta);
		list.push_back(nucleon);
		const ThreeVector &rcol = nucleon->getPosition();
		const ThreeVector zero;
		Particle *kaon = new Particle(KaonType,zero,rcol);
		list.push_back(kaon);
		
		PhaseSpaceGenerator::generateBiased(sqrtS, list, 0, angularSlope);
		
		
		fs->addModifiedParticle(delta);
		fs->addModifiedParticle(nucleon);
		fs->addCreatedParticle(kaon);
		
	}
}
