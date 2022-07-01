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

#include "G4INCLNpiToSKpiChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NpiToSKpiChannel::angularSlope = 6.; // What is the exact effect? Sould be check
	
	NpiToSKpiChannel::NpiToSKpiChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NpiToSKpiChannel::~NpiToSKpiChannel(){}
	
	void NpiToSKpiChannel::fillFinalState(FinalState *fs) {
		
        // p pi+ -> S+ pi+ K0 (5/4)
        // p pi+ -> S+ pi0 K+ (3/4)
        // p pi+ -> S0 pi+ K+ (1/4)
        //
        // p pi0 -> S+ pi0 K0 (1/2)
        // p pi0 -> S+ pi- K+ (1/2)
        // p pi0 -> S0 pi+ K0 (3/4)
        // p pi0 -> S0 pi0 K+ (3/8)
        // p pi0 -> S- pi+ K+ (1/2)
        //
        // p pi- -> S+ pi- K0 (3/8)
        // p pi- -> S0 pi0 K0 (5/8)
        // p pi- -> S0 pi- K+ (5/8)
        // p pi- -> S- pi+ K0 (1)
        // p pi- -> S- pi0 K+ (3/8)
		
		Particle *nucleon;
		Particle *pion;
		
		
		if(particle1->isNucleon()){
			nucleon = particle1;
			pion = particle2;
		}
		else{
			nucleon = particle2;
			pion = particle1;
		}
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(nucleon, pion);
		
		const G4int iso = ParticleTable::getIsospin(nucleon->getType()) + ParticleTable::getIsospin(pion->getType());
		G4double rdm = Random::shoot();
		
		ParticleType KaonType;
		
		if(iso == 3 || iso == -3){
			if(rdm*9. < 5.){
				KaonType = ParticleTable::getKaonType(-iso/3);
				nucleon->setType(ParticleTable::getSigmaType(iso/3*2));
			}
			else if(rdm*9. < 8.){
				KaonType = ParticleTable::getKaonType(iso/3);
				pion->setType(PiZero);
				nucleon->setType(ParticleTable::getSigmaType(iso/3*2));
			}
			else{
				KaonType = ParticleTable::getKaonType(iso/3);
				nucleon->setType(SigmaZero);
			}
		}
		else if(pion->getType() == PiZero){
			if(rdm*21. < 4.){
				KaonType = ParticleTable::getKaonType(-iso);
				nucleon->setType(ParticleTable::getSigmaType(iso*2));
			}
			else if(rdm*21. < 8.){
				KaonType = ParticleTable::getKaonType(iso);
				pion->setType(ParticleTable::getPionType(-2*iso));
				nucleon->setType(ParticleTable::getSigmaType(iso*2));
			}
			else if(rdm*21. < 14.){
				KaonType = ParticleTable::getKaonType(-iso);
				pion->setType(ParticleTable::getPionType(2*iso));
				nucleon->setType(SigmaZero);
			}
			else if(rdm*21. < 17.){
				KaonType = ParticleTable::getKaonType(iso);
				nucleon->setType(SigmaZero);
			}
			else{
				KaonType = ParticleTable::getKaonType(iso);
				pion->setType(ParticleTable::getPionType(2*iso));
				nucleon->setType(ParticleTable::getSigmaType(-iso*2));
			}
		}
		else{
			if(rdm*24. < 3.){
				KaonType = ParticleTable::getKaonType(iso);
				nucleon->setType(ParticleTable::getSigmaType(-iso*2));
			}
			else if(rdm*24. < 8.){
				KaonType = ParticleTable::getKaonType(iso);
				pion->setType(PiZero);
				nucleon->setType(SigmaZero);
			}
			else if(rdm*24. < 13.){
				KaonType = ParticleTable::getKaonType(-iso);
				nucleon->setType(SigmaZero);
			}
			else if(rdm*24. < 21.){
				KaonType = ParticleTable::getKaonType(iso);
				pion->setType(ParticleTable::getPionType(-2*iso));
				nucleon->setType(ParticleTable::getSigmaType(iso*2));
			}
			else{
				KaonType = ParticleTable::getKaonType(-iso);
				pion->setType(PiZero);
				nucleon->setType(ParticleTable::getSigmaType(iso*2));
			}
		}

		// Erase the parent resonance information of the nucleon and pion
		nucleon->setParentResonancePDGCode(0);
		nucleon->setParentResonanceID(0);
		pion->setParentResonancePDGCode(0);
		pion->setParentResonanceID(0);
		
		ParticleList list;
		list.push_back(nucleon);
		list.push_back(pion);
		const ThreeVector &rcol = nucleon->getPosition();
		const ThreeVector zero;
		Particle *kaon = new Particle(KaonType,zero,rcol);
		list.push_back(kaon);
		
		PhaseSpaceGenerator::generateBiased(sqrtS, list, 0, angularSlope);
		
		INCL_DEBUG("NpiToSKpi " << (kaon->getMomentum().theta()) * 180. / G4INCL::Math::pi << '\n');
		
		fs->addModifiedParticle(nucleon);
		fs->addModifiedParticle(pion);
		fs->addCreatedParticle(kaon);
		
	}
}
