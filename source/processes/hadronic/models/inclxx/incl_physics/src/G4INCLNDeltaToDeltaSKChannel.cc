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

#include "G4INCLNDeltaToDeltaSKChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NDeltaToDeltaSKChannel::angularSlope = 2.;
	
	NDeltaToDeltaSKChannel::NDeltaToDeltaSKChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NDeltaToDeltaSKChannel::~NDeltaToDeltaSKChannel(){}
	
	G4double NDeltaToDeltaSKChannel::sampleDeltaMass(G4double ecm) {
		const G4double maxDeltaMass = ecm - ParticleTable::effectiveSigmaMass - ParticleTable::effectiveKaonMass - 1.0;
		const G4double maxDeltaMassRndm = std::atan((maxDeltaMass-ParticleTable::effectiveDeltaMass)*2./ParticleTable::effectiveDeltaWidth); // atan((mass-1232)*2/130)
		const G4double deltaMassRndmRange = maxDeltaMassRndm - ParticleTable::minDeltaMassRndm; // atan
// assert(deltaMassRndmRange>0.);

		G4double y=ecm*ecm;
		G4double q2=(y-1.157776E6)*(y-6.4E5)/y/4.0; // 1.157776E6 = 1076^2 = (mNucleon + mPion)^2, 6.4E5 = 800^2 = (mNucleon - mPion)^2  // (prc56(1997)2431) (eq 3.7)^2
		G4double q3=std::pow(std::sqrt(q2), 3.);
		const G4double f3max=q3/(q3+5.832E6); // 5.832E6 = 180^3 = cut_parameter^3 // (prc56(1997)2431) (cf eq 3.6)
		G4double x;

		G4int nTries = 0;
		G4bool success = false;
		while(!success) { /* Loop checking, 10.07.2015, D.Mancusi */
			if(++nTries >= 100000) {
				INCL_WARN("NDeltaToDeltaSKChannel::sampleDeltaMass loop was stopped because maximum number of tries was reached. Minimum delta mass "
						  << ParticleTable::minDeltaMass << " MeV with CM energy " << ecm << " MeV may be unphysical." << '\n');
				return ParticleTable::minDeltaMass;
			}
			
			G4double rndm = ParticleTable::minDeltaMassRndm + Random::shoot() * deltaMassRndmRange; // atan in order to avec a distribution in 1/(1+x^2)
			y = std::tan(rndm); // (mass-1232)*2/130
			x = ParticleTable::effectiveDeltaMass + 0.5*ParticleTable::effectiveDeltaWidth*y; // probability to have the mass M = 1/(1+(M-1232)^2)/Pi cut with min and max mass
// assert(x>=ParticleTable::minDeltaMass && ecm >= x + ParticleTable::effectiveSigmaMass + ParticleTable::effectiveKaonMass + 1.0);
			
			// generation of the delta mass with the penetration factor
			// (see prc56(1997)2431)
			y=x*x;
			q2=(y-1.157776E6)*(y-6.4E5)/y/4.0; // 1.157776E6 = 1076^2 = (mNucleon + mPion)^2, 6.4E5 = 800^2 = (mNucleon - mPion)^2  // (prc56(1997)2431) (eq 3.7)^2
			q3=std::pow(std::sqrt(q2), 3.);
			const G4double f3=q3/(q3+5.832E6); // 5.832E6 = 180^3 = cut_parameter^3 // (prc56(1997)2431) (eq 3.6)
			rndm = Random::shoot();
			if (rndm*f3max < f3) success = true; // promoting high masses
		}
		return x;
	}
	
	void NDeltaToDeltaSKChannel::fillFinalState(FinalState *fs) {
        //
        // D++ p -> S+ K+ D+  (2)
        // D++ p -> S0 K+ D++ (1)
        // D++ p -> S+ K0 D++ (6)
        //
        // D++ n -> S+ K+ D0  (2)
        // D++ n -> S0 K+ D+  (4)
        // D++ n -> S- K+ D++ (6)
        // D++ n -> S+ K0 D+  (2)
        // D++ n -> S0 K0 D++ (1)
        //
        // D+  p -> S+ K+ D0  (2)
        // D+  p -> S0 K+ D+  (1)
        // D+  p -> S- K+ D++ (2)
        // D+  p -> S+ K0 D+  (2)
        // D+  p -> S0 K0 D++ (4)
        //
        // D+  n -> S+ K+ D-  (2)
        // D+  n -> S0 K+ D0  (4)
        // D+  n -> S- K+ D+  (2)
        // D+  n -> S+ K0 D0  (2)
        // D+  n -> S0 K0 D+  (1)
        // D+  n -> S- K0 D++ (2)
        
        
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
		ParticleType DeltaType;
		ParticleType SigmaType;
		
		const G4double rdm = Random::shoot();
		
		if(std::abs(iso) == 4){// D++ p
			if(rdm*9 < 2){
				KaonType = ParticleTable::getKaonType(iso/4);
				DeltaType = ParticleTable::getDeltaType(iso/4);
				SigmaType = ParticleTable::getSigmaType(iso/2);
			}
			else if(rdm*9 < 3){
				KaonType = ParticleTable::getKaonType(iso/4);
				DeltaType = ParticleTable::getDeltaType(3*iso/4);
				SigmaType = SigmaZero;
			}
			else{
				KaonType = ParticleTable::getKaonType(-iso/4);
				DeltaType = ParticleTable::getDeltaType(3*iso/4);
				SigmaType = ParticleTable::getSigmaType(iso/2);
			}
		}
		else if(iso == 0){// D+  n
			if(rdm*13 < 2){
				KaonType = ParticleTable::getKaonType(iso_d);
				DeltaType = ParticleTable::getDeltaType(-3*iso_d);
				SigmaType = ParticleTable::getSigmaType(2*iso_d);
			}
			else if(rdm*13 < 6){
				KaonType = ParticleTable::getKaonType(iso_d);
				DeltaType = ParticleTable::getDeltaType(-iso_d);
				SigmaType = SigmaZero;
			}
			else if(rdm*13 < 8){
				KaonType = ParticleTable::getKaonType(iso_d);
				DeltaType = ParticleTable::getDeltaType(iso_d);
				SigmaType = ParticleTable::getSigmaType(-2*iso_d);
			}
			else if(rdm*13 < 10){
				KaonType = ParticleTable::getKaonType(-iso_d);
				DeltaType = ParticleTable::getDeltaType(-iso_d);
				SigmaType = ParticleTable::getSigmaType(2*iso_d);
			}
			else if(rdm*13 < 11){
				KaonType = ParticleTable::getKaonType(-iso_d);
				DeltaType = ParticleTable::getDeltaType(iso_d);
				SigmaType = SigmaZero;
			}
			else{
				KaonType = ParticleTable::getKaonType(-iso_d);
				DeltaType = ParticleTable::getDeltaType(3*iso_d);
				SigmaType = ParticleTable::getSigmaType(-2*iso_d);
			}
		}
		else if(ParticleTable::getIsospin(particle1->getType()) == ParticleTable::getIsospin(particle2->getType())){// D+  p
			if(rdm*11 < 2){
				KaonType = ParticleTable::getKaonType(iso/2);
				DeltaType = ParticleTable::getDeltaType(-iso/2);
				SigmaType = ParticleTable::getSigmaType(iso);
			}
			else if(rdm*11 < 3){
				KaonType = ParticleTable::getKaonType(iso/2);
				DeltaType = ParticleTable::getDeltaType(iso/2);
				SigmaType = SigmaZero;
			}
			else if(rdm*11 < 5){
				KaonType = ParticleTable::getKaonType(iso/2);
				DeltaType = ParticleTable::getDeltaType(3*iso/2);
				SigmaType = ParticleTable::getSigmaType(-iso);
			}
			else if(rdm*11 < 7){
				KaonType = ParticleTable::getKaonType(-iso/2);
				DeltaType = ParticleTable::getDeltaType(iso/2);
				SigmaType = ParticleTable::getSigmaType(iso);
			}
			else{
				KaonType = ParticleTable::getKaonType(-iso/2);
				DeltaType = ParticleTable::getDeltaType(3*iso/2);
				SigmaType = SigmaZero;
			}
		}
		else{// D++ n 
			if(rdm*15 < 2){
				KaonType = ParticleTable::getKaonType(iso/2);
				DeltaType = ParticleTable::getDeltaType(-iso/2);
				SigmaType = ParticleTable::getSigmaType(iso);
			}
			else if(rdm*15 < 6){
				KaonType = ParticleTable::getKaonType(iso/2);
				DeltaType = ParticleTable::getDeltaType(iso/2);
				SigmaType = SigmaZero;
			}
			else if(rdm*15 < 12){
				KaonType = ParticleTable::getKaonType(iso/2);
				DeltaType = ParticleTable::getDeltaType(3*iso/2);
				SigmaType = ParticleTable::getSigmaType(-iso);
			}
			else if(rdm*15 < 14){
				KaonType = ParticleTable::getKaonType(-iso/2);
				DeltaType = ParticleTable::getDeltaType(iso/2);
				SigmaType = ParticleTable::getSigmaType(iso);
			}
			else{
				KaonType = ParticleTable::getKaonType(-iso/2);
				DeltaType = ParticleTable::getDeltaType(3*iso/2);
				SigmaType = SigmaZero;
			}
		}
		
				
		particle1->setType(DeltaType);
		delta->setMass(sampleDeltaMass(sqrtS));
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
