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

#include "G4INCLNNToNLKpiChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NNToNLKpiChannel::angularSlope = 2.; // What is the exact effect? Sould be check
	
	NNToNLKpiChannel::NNToNLKpiChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NNToNLKpiChannel::~NNToNLKpiChannel(){}
	
	void NNToNLKpiChannel::fillFinalState(FinalState *fs) {
		
        // pp -> p pi+ L K0 (9)
        // pp -> p pi0 L K+ (2)
        // pp -> n pi+ L K+ (1)
        // pn -> p pi- L K+ (2)
        // pn -> n pi0 L K+ (4)
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(particle1, particle2);
		
		const G4int iso = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
		
		ParticleType KaonType;
		ParticleType PionType;
		
		G4double rdm = Random::shoot();
		particle2->setType(Lambda);
		
		if(iso == 2){
			if(rdm*4. < 3.){
				KaonType = KZero;
				PionType = PiPlus;
			}
			else if(rdm*12. < 11.){
				KaonType = KPlus;
				PionType = PiZero;
			}
			else{
				KaonType = KPlus;
				PionType = PiPlus;
				particle1->setType(Neutron);
			}
		}
		else if(iso == -2){
			if(rdm*4. < 3.){
				KaonType = KPlus;
				PionType = PiMinus;
			}
			else if(rdm*12. < 11.){
				KaonType = KZero;
				PionType = PiZero;
			}
			else{
				KaonType = KZero;
				PionType = PiMinus;
				particle1->setType(Proton);
			}
		}
		else if(rdm*6. < 2.){
			KaonType = KPlus;
			PionType = PiZero;
			particle1->setType(Neutron);
		}
		else if(rdm*6. < 4.){
			KaonType = KZero;
			PionType = PiZero;
			particle1->setType(Proton);
		}
		else if(rdm*6. < 5.){
			KaonType = KPlus;
			PionType = PiMinus;
			particle1->setType(Proton);
		}
		else{
			KaonType = KZero;
			PionType = PiPlus;
			particle1->setType(Neutron);
		}
		
		ParticleList list;
		list.push_back(particle1);
		list.push_back(particle2);
		const ThreeVector &rcol1 = particle1->getPosition();
		const ThreeVector &rcol2 = particle2->getPosition();
		const ThreeVector zero;
		Particle *pion = new Particle(PionType,zero,rcol1);
		Particle *kaon = new Particle(KaonType,zero,rcol2);
		list.push_back(kaon);
		list.push_back(pion);
		
		PhaseSpaceGenerator::generateBiased(sqrtS, list, 0, angularSlope);
		
		INCL_DEBUG("NNToNLKpi " << (kaon->getMomentum().theta()) * 180. / G4INCL::Math::pi << '\n');
		
		fs->addModifiedParticle(particle1);
		fs->addModifiedParticle(particle2);
		fs->addCreatedParticle(kaon);
		fs->addCreatedParticle(pion);
		
	}
}
