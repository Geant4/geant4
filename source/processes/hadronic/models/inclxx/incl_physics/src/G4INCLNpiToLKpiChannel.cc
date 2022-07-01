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

#include "G4INCLNpiToLKpiChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NpiToLKpiChannel::angularSlope = 6.; // What is the exact effect? Sould be check
	
	NpiToLKpiChannel::NpiToLKpiChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NpiToLKpiChannel::~NpiToLKpiChannel(){}
	
	void NpiToLKpiChannel::fillFinalState(FinalState *fs) {
		
        // ratio
        // p pi+ (1)	p pi0 (3/2)		p pi- (2)
        //
        // p pi0 -> L K+ pi0 (1/2)
        // all the others (1)
        
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
			KaonType = ParticleTable::getKaonType(iso/3);
		}
		else if(pion->getType() == PiZero){
			if(rdm*3. < 1.){
				KaonType = ParticleTable::getKaonType(iso);
			}
			else{
				KaonType = ParticleTable::getKaonType(-iso);
				pion->setType(ParticleTable::getPionType(2*iso));
			}
		}
		else{
			if(rdm < 0.5){
				KaonType = ParticleTable::getKaonType(-iso);
			}
			else{
				KaonType = ParticleTable::getKaonType(iso);
				pion->setType(PiZero);
			}
		}
		
		nucleon->setType(Lambda);

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
		
		INCL_DEBUG("NpiToLKpi " << (kaon->getMomentum().theta()) * 180. / G4INCL::Math::pi << '\n');
		
		fs->addModifiedParticle(nucleon);
		fs->addModifiedParticle(pion);
		fs->addCreatedParticle(kaon);
		
	}
}
