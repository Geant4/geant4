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

#include "G4INCLNKToNKpiChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NKToNKpiChannel::angularSlope = 4.; // What is the exact effect? Sould be check
	
	NKToNKpiChannel::NKToNKpiChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NKToNKpiChannel::~NKToNKpiChannel(){}
	
	void NKToNKpiChannel::fillFinalState(FinalState *fs) {
		
        // K+ p -> p K+ pi0		(1.2)
        // K+ p -> p K0 pi+		(3)
        // K+ p -> n K+ pi+		(0.8)
        //
        // K0 p -> p K+ pi-		(1)
        // K0 p -> p K0 pi0		(0.845)
        // K0 p -> n K+ pi0		(1.47)
        // K0 p -> n K0 pi+		(2.23)
        
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
			if(rdm*5. < 1.2){
				PionType = PiZero;
			}
			else if(rdm*5. < 4.2){
				PionType = ParticleTable::getPionType(iso);
				kaon->setType(ParticleTable::getKaonType(-iso/2));
			}
			else{
				PionType = ParticleTable::getPionType(iso);
				nucleon->setType(ParticleTable::getNucleonType(-iso/2));
			}
		}
		else{
			if(rdm*5.545 < 1.){
				kaon->setType(ParticleTable::getKaonType(iso_n));
				PionType = ParticleTable::getPionType(-2*iso_n);
			}
			else if(rdm*5.545 < 1.845){
				PionType = PiZero;
			}
			else if(rdm*5.545 < 3.315){
				kaon->setType(ParticleTable::getKaonType(iso_n));
				nucleon->setType(ParticleTable::getNucleonType(-iso_n));
				PionType = PiZero;
			}
			else{
				nucleon->setType(ParticleTable::getNucleonType(-iso_n));
				PionType = ParticleTable::getPionType(2*iso_n);
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
