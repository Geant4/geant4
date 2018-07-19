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

#include "G4INCLNNToNLK2piChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NNToNLK2piChannel::angularSlope = 2.; // What is the exact effect? Sould be check
	
	NNToNLK2piChannel::NNToNLK2piChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NNToNLK2piChannel::~NNToNLK2piChannel(){}
	
	void NNToNLK2piChannel::fillFinalState(FinalState *fs) {
		
		/* Equipartition in all channel with factor N(pi)!
		*/
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(particle1, particle2);
		
		const G4int iso = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
		
		ParticleType KaonType;
		ParticleType Pion1Type;
		ParticleType Pion2Type;
		
		G4double rdm = Random::shoot();
		particle2->setType(Lambda);
		
		if(iso == 2){
			if(rdm*7. < 2.){
				particle1->setType(Neutron);
				KaonType =  KZero;
				Pion1Type =  PiPlus;
				Pion2Type =  PiPlus;
			}
			else if(rdm*7. < 3.){
				particle1->setType(Neutron);
				KaonType =  KPlus;
				Pion1Type =  PiZero;
				Pion2Type =  PiPlus;
			}
			else if(rdm*7. < 4.){
				particle1->setType(Proton);
				KaonType =  KZero;
				Pion1Type =  PiZero;
				Pion2Type =  PiPlus;
			}
			else if(rdm*7. < 5.){
				particle1->setType(Proton);
				KaonType =  KPlus;
				Pion1Type =  PiMinus;
				Pion2Type =  PiPlus;
			}
			else{
				particle1->setType(Proton);
				KaonType =  KPlus;
				Pion1Type =  PiZero;
				Pion2Type =  PiZero;
			}
			
		}if(iso == -2){
			if(rdm*7. < 1.){
				particle1->setType(Neutron);
				KaonType =  KZero;
				Pion1Type =  PiMinus;
				Pion2Type =  PiPlus;
			}
			else if(rdm*7. < 3.){
				particle1->setType(Neutron);
				KaonType =  KZero;
				Pion1Type =  PiZero;
				Pion2Type =  PiZero;
			}
			else if(rdm*7. < 4.){
				particle1->setType(Neutron);
				KaonType =  KPlus;
				Pion1Type =  PiMinus;
				Pion2Type =  PiZero;
			}
			else if(rdm*7. < 5.){
				particle1->setType(Proton);
				KaonType =  KZero;
				Pion1Type =  PiMinus;
				Pion2Type =  PiZero;
			}
			else{
				particle1->setType(Proton);
				KaonType =  KPlus;
				Pion1Type =  PiMinus;
				Pion2Type =  PiMinus;
			}
		}
		else{
			if(rdm*8. < 1.){
				particle1->setType(Neutron);
				KaonType =  KZero;
				Pion1Type =  PiZero;
				Pion2Type =  PiPlus;
			}
			else if(rdm*8. < 2.){
				particle1->setType(Neutron);
				KaonType =  KPlus;
				Pion1Type =  PiMinus;
				Pion2Type =  PiPlus;
			}
			else if(rdm*8. < 4.){
				particle1->setType(Neutron);
				KaonType =  KPlus;
				Pion1Type =  PiZero;
				Pion2Type =  PiZero;
			}
			else if(rdm*8. < 5.){
				particle1->setType(Proton);
				KaonType =  KZero;
				Pion1Type =  PiMinus;
				Pion2Type =  PiPlus;
			}
			else if(rdm*8. < 7.){
				particle1->setType(Proton);
				KaonType =  KZero;
				Pion1Type =  PiZero;
				Pion2Type =  PiZero;
			}
			else{
				particle1->setType(Proton);
				KaonType =  KPlus;
				Pion1Type =  PiMinus;
				Pion2Type =  PiZero;
			}
		}
		
		
		ParticleList list;
		list.push_back(particle1);
		list.push_back(particle2);
		const ThreeVector &rcol1 = particle1->getPosition();
		const ThreeVector &rcol2 = particle2->getPosition();
		const ThreeVector zero;
		Particle *pion1 = new Particle(Pion1Type,zero,rcol1);
		Particle *pion2 = new Particle(Pion2Type,zero,rcol1);
		Particle *kaon = new Particle(KaonType,zero,rcol2);
		list.push_back(kaon);
		list.push_back(pion1);
		list.push_back(pion2);
		
		PhaseSpaceGenerator::generateBiased(sqrtS, list, 0, angularSlope);
		
		INCL_DEBUG("NNToNLK2pi " << (kaon->getMomentum().theta()) * 180. / G4INCL::Math::pi << '\n');
		
		fs->addModifiedParticle(particle1);
		fs->addModifiedParticle(particle2);
		fs->addCreatedParticle(kaon);
		fs->addCreatedParticle(pion1);
		fs->addCreatedParticle(pion2);
		
	}
}
