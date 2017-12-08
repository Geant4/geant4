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

#include "G4INCLNKbToNKb2piChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NKbToNKb2piChannel::angularSlope = 4.; // What is the exact effect? Sould be check
	
	NKbToNKb2piChannel::NKbToNKb2piChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NKbToNKb2piChannel::~NKbToNKb2piChannel(){}
	
	void NKbToNKb2piChannel::fillFinalState(FinalState *fs) {
		
        // p K0b -> p K0b pi+ pi- (1)
        // p K0b -> p K0b pi0 pi0 (1/4)
        // p K0b -> p K-  pi+ pi0 (1)
        // p K0b -> n K0b pi+ pi0 (1)
        // p K0b -> n K-  pi+ pi+ (1)
        //
        // p K-  -> p K-  pi+ pi- (1)
        // p K-  -> p K-  pi0 pi0 (1/4)
        // p K-  -> p K0b pi0 pi- (1)
        // p K-  -> n K0b pi+ pi- (1)
        // p K-  -> n K0b pi0 pi0 (1/2)
        // p K-  -> n K-  pi+ pi0 (1)
        
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
		
		ParticleType Pion1Type;
		ParticleType Pion2Type;
		
		if(iso == 2 || iso == -2){
			if(rdm*17. < 4.){
				Pion1Type = PiPlus;
				Pion2Type = PiMinus;
			}
			else if(rdm*17. < 5.){
				Pion1Type = PiZero;
				Pion2Type = PiZero;
			}
			else if(rdm*17. < 9.){
				Pion1Type = ParticleTable::getPionType(iso);
				Pion2Type = PiZero;
				kaon->setType(ParticleTable::getAntiKaonType(-iso_n));
			}
			else if(rdm*17. < 13.){
				Pion1Type = ParticleTable::getPionType(iso);
				Pion2Type = PiZero;
				nucleon->setType(ParticleTable::getNucleonType(-iso_n));
			}
			else{
				Pion1Type = ParticleTable::getPionType(iso);
				Pion2Type = ParticleTable::getPionType(iso);
				kaon->setType(ParticleTable::getAntiKaonType(-iso_n));
				nucleon->setType(ParticleTable::getNucleonType(-iso_n));
			}
		}
		else{
			if(rdm*19. < 4.){
				Pion1Type = PiPlus;
				Pion2Type = PiMinus;
			}
			else if(rdm*19. < 5.){
				Pion1Type = PiZero;
				Pion2Type = PiZero;
			}
			else if(rdm*19. < 9.){
				Pion1Type = ParticleTable::getPionType(-2*iso_n);
				Pion2Type = PiZero;
				kaon->setType(ParticleTable::getAntiKaonType(iso_n));
			}
			else if(rdm*19. < 13.){
				Pion1Type = PiPlus;
				Pion2Type = PiMinus;
				kaon->setType(ParticleTable::getAntiKaonType(iso_n));
				nucleon->setType(ParticleTable::getNucleonType(-iso_n));
			}
			else if(rdm*19. < 15.){
				Pion1Type = PiZero;
				Pion2Type = PiZero;
				kaon->setType(ParticleTable::getAntiKaonType(iso_n));
				nucleon->setType(ParticleTable::getNucleonType(-iso_n));
			}
			else{
				Pion1Type = ParticleTable::getPionType(2*iso_n);
				Pion2Type = PiZero;
				nucleon->setType(ParticleTable::getNucleonType(-iso_n));
			}
		}
		
		ParticleList list;
		list.push_back(nucleon);
		list.push_back(kaon);
		const ThreeVector &rcol1 = nucleon->getPosition();
		const ThreeVector &rcol2 = kaon->getPosition();
		const ThreeVector zero;
		Particle *pion1 = new Particle(Pion1Type,zero,rcol1);
		Particle *pion2 = new Particle(Pion2Type,zero,rcol2);
		list.push_back(pion1);
		list.push_back(pion2);
		
		PhaseSpaceGenerator::generateBiased(sqrtS, list, 0, angularSlope);
		
		fs->addModifiedParticle(nucleon);
		fs->addModifiedParticle(kaon);
		fs->addCreatedParticle(pion1);
		fs->addCreatedParticle(pion2);
		
	}
}
