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

#include "G4INCLNKbToS2piChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NKbToS2piChannel::angularSlope = 4.; // What is the exact effect? Sould be check
	
	NKbToS2piChannel::NKbToS2piChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NKbToS2piChannel::~NKbToS2piChannel(){}
	
	void NKbToS2piChannel::fillFinalState(FinalState *fs) {
        
        // p K0b -> S+ pi+ pi- (2/3)
        // p K0b -> S+ pi0 pi0 (1/4)
        // p K0b -> S0 pi+ pi0 (5/6)
        // p K0b -> S- pi+ pi+ (2/3)
        //
        // p K-  -> S+ pi0 pi- (1)
        // p K-  -> S0 pi+ pi- (2/3)
        // p K-  -> S0 pi0 pi0 (1/8)
        // p K-  -> S- pi+ pi0 (2/3)
        
		Particle *nucleon;
		Particle *kaon;
				
		if(particle1->isNucleon()){
			nucleon = particle1;
			kaon = particle2;
		}
		else{
			nucleon = particle2;
			kaon = particle1;
		}
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(nucleon, kaon);
		
		const G4int iso = ParticleTable::getIsospin(nucleon->getType()) + ParticleTable::getIsospin(kaon->getType());
		const G4int iso_n = ParticleTable::getIsospin(nucleon->getType());
		G4double rdm = Random::shoot();
		
		ParticleType PionType;
		
		if(iso == 2 || iso == -2){
			if(rdm*29. < 8.){
				PionType = ParticleTable::getPionType(-iso);
				kaon->setType(ParticleTable::getPionType(iso));
				nucleon->setType(ParticleTable::getSigmaType(iso));
			}
			else if(rdm*29. < 11.){
				PionType = PiZero;
				kaon->setType(PiZero);
				nucleon->setType(ParticleTable::getSigmaType(iso));
			}
			else if(rdm*29. < 21.){
				PionType = PiZero;
				kaon->setType(ParticleTable::getPionType(iso));
				nucleon->setType(SigmaZero);
			}
			else{
				PionType = ParticleTable::getPionType(iso);
				kaon->setType(ParticleTable::getPionType(iso));
				nucleon->setType(ParticleTable::getSigmaType(-iso));
			}
		}
		else{
			if(rdm*59. < 24.){
				PionType = PiZero;
				kaon->setType(ParticleTable::getPionType(-2*iso_n));
				nucleon->setType(ParticleTable::getSigmaType(2*iso_n));
			}
			else if(rdm*59. < 40.){
				PionType = ParticleTable::getPionType(2*iso_n);
				kaon->setType(ParticleTable::getPionType(-2*iso_n));
				nucleon->setType(SigmaZero);
			}
			else if(rdm*59. < 43.){
				PionType = PiZero;
				kaon->setType(PiZero);
				nucleon->setType(SigmaZero);
			}
			else{
				PionType = ParticleTable::getPionType(2*iso_n);
				kaon->setType(PiZero);
				nucleon->setType(ParticleTable::getSigmaType(-2*iso_n));
			}
		}
		
		ParticleList list;
		list.push_back(nucleon);
		list.push_back(kaon);
		const ThreeVector &rcol = nucleon->getPosition();
		const ThreeVector zero;
		Particle *pion = new Particle(PionType,zero,rcol);
		list.push_back(pion);
		
		PhaseSpaceGenerator::generateBiased(sqrtS, list, 0, angularSlope);
		
		fs->addModifiedParticle(nucleon);
		fs->addModifiedParticle(kaon);
		fs->addCreatedParticle(pion);
		
	}
}
