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

#include "G4INCLNpiToMissingStrangenessChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NpiToMissingStrangenessChannel::angularSlope = 1.;
	
	NpiToMissingStrangenessChannel::NpiToMissingStrangenessChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NpiToMissingStrangenessChannel::~NpiToMissingStrangenessChannel(){}
	
	void NpiToMissingStrangenessChannel::fillFinalState(FinalState *fs) {
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(particle1, particle2); // MeV		/!\/!\/!\.
// assert(sqrtS > 2.240); // ! > 2.109 Not supposed to be under 2.244 GeV.
		
		const G4int iso = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
// assert(iso == -3 || iso == -1 || iso == 1 || iso == 3);
		G4int iso_system = 0;
		G4int available_iso = 0;
		G4int nbr_pions = 0;
		G4int min_pions = 0;
		G4int max_pions = 0;
		
		Particle *pion_initial;
		Particle *nucleon_initial;
		
		if(particle1->isPion()){
			pion_initial = particle1;
			nucleon_initial = particle2;
		}
		else{
			pion_initial = particle2;
			nucleon_initial = particle1;
		}
		const G4double pLab = 0.001*KinematicsUtils::momentumInLab(pion_initial, nucleon_initial); // GeV	/!\/!\/!\.
		
		G4double rdm = Random::shoot();
		
		if(rdm < 0.35){
			// Lambda-K chosen
			nucleon_initial->setType(Lambda);
			available_iso = 1;
			min_pions = 3;
			max_pions = G4int((sqrtS-ParticleTable::getINCLMass(Lambda)-ParticleTable::getINCLMass(KZero)-10.)/ParticleTable::getINCLMass(PiPlus));
		}
		else if((iso == 0 && rdm < 0.55) || rdm < 0.5){
			// N-K-Kb chosen
			available_iso = 3;
			min_pions = 1;
			max_pions = G4int((sqrtS-ParticleTable::getINCLMass(Proton)-2.*ParticleTable::getINCLMass(KZero)-10.)/ParticleTable::getINCLMass(PiPlus));
		}
		else{
			// Sigma-K chosen
			available_iso = 3;
			min_pions = 3;
			max_pions = G4int((sqrtS-ParticleTable::getINCLMass(SigmaMinus)-ParticleTable::getINCLMass(KZero)-10.)/ParticleTable::getINCLMass(PiPlus));
		}
										// Gaussian noise  + mean value nbr pions fonction energy (choice)
		G4double intermediaire = min_pions + Random::gauss(2) + std::sqrt(pLab-2.2);
		nbr_pions = std::min(max_pions,std::max(min_pions,G4int(intermediaire )));
		
		available_iso += nbr_pions*2;

		// Erase the parent resonance information of the initial particles
		particle1->setParentResonancePDGCode(0);
		particle1->setParentResonanceID(0);
		particle2->setParentResonancePDGCode(0);
		particle2->setParentResonanceID(0);
		
		ParticleList list;
		ParticleType PionType = PiZero;
		const ThreeVector &rcol1 = pion_initial->getPosition();
		const ThreeVector zero;
		
		//      (pip   piz   pim)   (sp    sz    sm)   (L    S    Kb)
		//pip_p  0.63  0.26  0.11   0.73  0.25  0.02   0.42 0.49 0.09 // inital
		//pip_p  0.54  0.26  0.20   0.73  0.25  0.02   0.42 0.49 0.09 // choice
		G4bool pip_p = (std::abs(iso) == 3);
		//piz_p  0.32  0.45  0.23   0.52  0.40  0.08   0.40 0.41 0.19
		G4bool piz_p = (ParticleTable::getIsospin(pion_initial->getType()) == 0);
		//pim_p  0.18  0.37  0.45   0.20  0.63  0.17   0.39 0.33 0.28
		G4bool pim_p = (!pip_p && !piz_p);
	
		for(Int_t i=0; i<nbr_pions; i++){
			Particle *pion = new Particle(PionType,zero,rcol1);
			if(available_iso-std::abs(iso-iso_system) >= 4){
				rdm = Random::shoot();
				if((pip_p && rdm < 0.54) || (piz_p && rdm < 0.32) || (pim_p && rdm < 0.45)){
					pion->setType(ParticleTable::getPionType(G4int(Math::sign(iso))*2)); //pip/pip/pim
					iso_system += 2*G4int(Math::sign(iso));
					available_iso -= 2;
				}
				else if((pip_p && rdm < 0.80) || (piz_p && rdm < 0.77) || (pim_p && rdm < 0.82)){
					pion->setType(PiZero);
					available_iso -= 2;
				}
				else{
					pion->setType(ParticleTable::getPionType(-G4int(Math::sign(iso))*2));
					iso_system -= 2*G4int(Math::sign(iso));
					available_iso -= 2;
				}
			}
			else if(available_iso-std::abs(iso-iso_system) == 2){
				rdm = Random::shoot();
				if((pip_p && rdm < 0.26/0.37 && (Math::sign(iso)*Math::sign(iso-iso_system)+1)) || (pip_p && rdm < 0.26/0.89 && (Math::sign(iso)*Math::sign(iso-iso_system)-1)) ||
				   (piz_p && rdm < 0.45/0.68 && (Math::sign(iso)*Math::sign(iso-iso_system)+1)) || (piz_p && rdm < 0.45/0.77 && (Math::sign(iso)*Math::sign(iso-iso_system)-1)) ||
				   (pim_p && rdm < 0.37/0.82 && (Math::sign(iso)*Math::sign(iso-iso_system)+1)) || (piz_p && rdm < 0.37/0.55 && (Math::sign(iso)*Math::sign(iso-iso_system)-1))){
					pion->setType(PiZero);
					available_iso -= 2;
				}
				else{
					pion->setType(ParticleTable::getPionType(Math::sign(iso-iso_system)*2));
					iso_system += Math::sign(iso-iso_system)*2;
					available_iso -= 2;
				}
			}
			else if(available_iso-std::abs(iso-iso_system) == 0){
				pion->setType(ParticleTable::getPionType(Math::sign(iso-iso_system)*2));
				iso_system += Math::sign(iso-iso_system)*2;
				available_iso -= 2;
			}
			else INCL_ERROR("Pion Generation Problem in NpiToMissingStrangenessChannel" << '\n');
			list.push_back(pion);
		}
		
		if(nucleon_initial->isLambda()){ // K-Lambda
// assert(available_iso == 1);
			pion_initial->setType(ParticleTable::getKaonType(iso-iso_system));
		}
		else if(min_pions == 1){ // N-K-Kb chosen
// assert(available_iso == 3);
			Particle *antikaon = new Particle(KMinus,zero,rcol1);
			if(std::abs(iso-iso_system) == 3){
				pion_initial->setType(ParticleTable::getKaonType((iso-iso_system)/3));
				nucleon_initial->setType(ParticleTable::getNucleonType((iso-iso_system)/3));
				antikaon->setType(ParticleTable::getAntiKaonType((iso-iso_system)/3));
			}
			else if(std::abs(iso-iso_system) == 1){ // equi-repartition
				rdm = G4int(Random::shoot()*3.)-1;
				nucleon_initial->setType(ParticleTable::getNucleonType((G4int(rdm+0.5)*2-1)*(iso_system-iso)));
				pion_initial->setType(ParticleTable::getKaonType((std::abs(rdm*2)-1)*(iso-iso_system)));
				antikaon->setType(ParticleTable::getAntiKaonType((G4int(rdm-0.5)*2+1)*(iso-iso_system)));
			}
			else INCL_ERROR("Isospin non-conservation in NNToMissingStrangenessChannel" << '\n');
			list.push_back(antikaon);
			nbr_pions += 1; // need for addCreatedParticle loop
		}
		else{// Sigma-K
// assert(available_iso == 3);
			if(std::abs(iso-iso_system) == 3){
				pion_initial->setType(ParticleTable::getKaonType((iso-iso_system)/3));
				nucleon_initial->setType(ParticleTable::getSigmaType((iso-iso_system)*2/3));
			}
			else if(std::abs(iso-iso_system) == 1){
				rdm = Random::shoot();
				if((pip_p && rdm < 0.73) || (piz_p && rdm < 0.32) || (pim_p && rdm < 0.45)){
					nucleon_initial->setType(SigmaZero);
					pion_initial->setType(ParticleTable::getKaonType(iso-iso_system));
				}
				else{
					nucleon_initial->setType(ParticleTable::getSigmaType((iso-iso_system)*2));
					pion_initial->setType(ParticleTable::getKaonType(iso_system-iso));
				}
			}
			else INCL_ERROR("Isospin non-conservation in NNToMissingStrangenessChannel" << '\n');
		}
		
		list.push_back(pion_initial);
		list.push_back(nucleon_initial);
		
		PhaseSpaceGenerator::generateBiased(sqrtS, list, list.size()-1, angularSlope);
		
		fs->addModifiedParticle(pion_initial);
		fs->addModifiedParticle(nucleon_initial);
		for(Int_t i=0; i<nbr_pions; i++) fs->addCreatedParticle(list[i]);
		
	}
}
