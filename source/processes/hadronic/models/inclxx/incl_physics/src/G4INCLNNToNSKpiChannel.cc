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

#include "G4INCLNNToNSKpiChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NNToNSKpiChannel::angularSlope = 2.; // What is the exact effect? Sould be check
	
	NNToNSKpiChannel::NNToNSKpiChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NNToNSKpiChannel::~NNToNSKpiChannel(){}
	
	void NNToNSKpiChannel::fillFinalState(FinalState *fs) {
		
		// pp (36)	pn (36)
		//
        // pp -> p pi+ S- K+ (9)
        // pp -> p pi+ S0 K0 (9)
        // pp -> p pi0 S+ K0 (4)
        // pp -> n pi+ S+ K0 (2)
        // pp -> p pi0 S0 K+ (4)
        // pp -> n pi+ S0 K+ (2)
        // pp -> p pi- S+ K+ (2)
        // pp -> n pi0 S+ K+ (4)
        
        // pn -> p pi0 S- K+ (4)
        // pn -> n pi+ S- K+ (2)
        // pn -> p pi0 S0 K0 (2)
        // pn -> n pi+ S0 K0 (1)
        // pn -> p pi+ S- K0 (9)
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(particle1, particle2);
		
		const G4int iso = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
		
		ParticleType KaonType;
		ParticleType PionType;
		
		G4double rdm = Random::shoot();
		
		if(iso == 2){
			if(rdm * 36. < 9.){
				KaonType = KPlus;
				PionType = PiPlus;
				particle2->setType(SigmaMinus);
			}
			else if(rdm * 36. < 18.){
				KaonType = KZero;
				PionType = PiPlus;
				particle2->setType(SigmaZero);
			}
			else if(rdm * 36. < 22.){
				KaonType = KZero;
				PionType = PiZero;
				particle2->setType(SigmaPlus);
			}
			else if(rdm * 36. < 24.){
				KaonType = KZero;
				PionType = PiPlus;
				particle1->setType(Neutron);
				particle2->setType(SigmaPlus);
			}
			else if(rdm * 36. < 28.){
				KaonType = KPlus;
				PionType = PiZero;
				particle2->setType(SigmaZero);
			}
			else if(rdm * 36. < 30.){
				KaonType = KPlus;
				PionType = PiPlus;
				particle1->setType(Neutron);
				particle2->setType(SigmaZero);
			}
			else if(rdm * 36. < 32.){
				KaonType = KPlus;
				PionType = PiMinus;
				particle2->setType(SigmaPlus);
			}
			else{
				KaonType = KPlus;
				PionType = PiZero;
				particle1->setType(Neutron);
				particle2->setType(SigmaPlus);
			}
			
		}
		else if(iso == -2){
			if(rdm * 36. < 9.){
				KaonType = KZero;
				PionType = PiMinus;
				particle2->setType(SigmaPlus);
			}
			else if(rdm * 36. < 18.){
				KaonType = KPlus;
				PionType = PiMinus;
				particle2->setType(SigmaZero);
			}
			else if(rdm * 36. < 22.){
				KaonType = KPlus;
				PionType = PiZero;
				particle2->setType(SigmaMinus);
			}
			else if(rdm * 36. < 24.){
				KaonType = KPlus;
				PionType = PiMinus;
				particle1->setType(Proton);
				particle2->setType(SigmaMinus);
			}
			else if(rdm * 36. < 28.){
				KaonType = KZero;
				PionType = PiZero;
				particle2->setType(SigmaZero);
			}
			else if(rdm * 36. < 30.){
				KaonType = KZero;
				PionType = PiMinus;
				particle1->setType(Proton);
				particle2->setType(SigmaZero);
			}
			else if(rdm * 36. < 32.){
				KaonType = KZero;
				PionType = PiPlus;
				particle2->setType(SigmaMinus);
			}
			else{
				KaonType = KZero;
				PionType = PiZero;
				particle1->setType(Proton);
				particle2->setType(SigmaMinus);
			}
			
		}
		else if(rdm*36. < 4.){
			KaonType = KPlus;
			PionType = PiZero;
			particle1->setType(Proton);
			particle2->setType(SigmaMinus);
		}
		else if(rdm*36. < 6.){
			KaonType = KZero;
			PionType = PiZero;
			particle1->setType(Neutron);
			particle2->setType(SigmaPlus);
		}
		else if(rdm*36. < 8.){
			KaonType = KPlus;
			PionType = PiPlus;
			particle1->setType(Neutron);
			particle2->setType(SigmaMinus);
		}
		else if(rdm*36. < 9.){
			KaonType = KZero;
			PionType = PiMinus;
			particle1->setType(Proton);
			particle2->setType(SigmaPlus);
		}
		else if(rdm*36. < 18.){
			KaonType = KZero;
			PionType = PiZero;
			particle1->setType(Proton);
			particle2->setType(SigmaZero);
		}
		else if(rdm*36. < 27.){
			KaonType = KPlus;
			PionType = PiZero;
			particle1->setType(Neutron);
			particle2->setType(SigmaZero);
		}
		else if(rdm*36. < 28.){
			KaonType = KZero;
			PionType = PiPlus;
			particle1->setType(Neutron);
			particle2->setType(SigmaZero);
		}
		else if(rdm*36. < 30.){
			KaonType = KPlus;
			PionType = PiMinus;
			particle1->setType(Proton);
			particle2->setType(SigmaZero);
		}
		else if(rdm*36. < 32.){
			KaonType = KZero;
			PionType = PiPlus;
			particle1->setType(Proton);
			particle2->setType(SigmaMinus);
		}
		else{
			KaonType = KPlus;
			PionType = PiMinus;
			particle1->setType(Neutron);
			particle2->setType(SigmaPlus);
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
		
		INCL_DEBUG("NNToNSKpi " << (kaon->getMomentum().theta()) * 180. / G4INCL::Math::pi << '\n');
		
		fs->addModifiedParticle(particle1);
		fs->addModifiedParticle(particle2);
		fs->addCreatedParticle(kaon);
		fs->addCreatedParticle(pion);
		
	}
}
