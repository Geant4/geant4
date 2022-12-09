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

#include "G4INCLNNToMissingStrangenessChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NNToMissingStrangenessChannel::angularSlope = 1.;
	
	NNToMissingStrangenessChannel::NNToMissingStrangenessChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NNToMissingStrangenessChannel::~NNToMissingStrangenessChannel(){}
	
	void NNToMissingStrangenessChannel::fillFinalState(FinalState *fs) {
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(particle1, particle2); // MeV
		const G4double pLab = 0.001*KinematicsUtils::momentumInLab(particle1, particle2); // GeV
// assert(sqrtS > 3100); // ! > 3047.34. Not supposed to be under 3,626 GeV.
		
		const G4int iso = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
// assert(iso == -2 || iso == 0 || iso == 2);
		G4int iso_system = 0;
		G4int available_iso = 0;
		G4int nbr_pions = 0;
		G4int min_pions = 0;
		G4int max_pions = 0;
		
		G4double rdm = Random::shoot();
		
		if(rdm < 0.35){
			// N-K-Lambda chosen
			particle2->setType(Lambda);
			available_iso = 2;
			min_pions = 3;
			max_pions = G4int((sqrtS-ParticleTable::getINCLMass(Lambda)-ParticleTable::getINCLMass(Proton)-ParticleTable::getINCLMass(KZero))/ParticleTable::getINCLMass(PiPlus));
		}
		else if((iso == 0 && rdm < 0.55) || rdm < 0.5){
			// N-N-K-Kb chosen
			available_iso = 4;
			min_pions = 1;
			max_pions = G4int((sqrtS-2.*ParticleTable::getINCLMass(Proton)-2.*ParticleTable::getINCLMass(KZero))/ParticleTable::getINCLMass(PiPlus));
		}
		else{
			// N-K-Sigma chosen
			available_iso = 4;
			min_pions = 3;
			max_pions = G4int((sqrtS-ParticleTable::getINCLMass(SigmaMinus)-ParticleTable::getINCLMass(Proton)-ParticleTable::getINCLMass(KZero))/ParticleTable::getINCLMass(PiPlus));
		}
										/* Gaussian noise  + mean value nbr pions fonction energy (choice)*/
		G4double intermediaire = min_pions + Random::gauss(2) + std::sqrt(pLab-5);
		nbr_pions = std::min(max_pions,std::max(min_pions,G4int(intermediaire )));
		
		available_iso += nbr_pions*2;
		
		ParticleList list;
		ParticleType PionType = PiZero;
		const ThreeVector &rcol1 = particle1->getPosition();
		const ThreeVector zero;
		Particle *kaon = new Particle(KPlus,zero,rcol1);
		
		for(Int_t i=0; i<nbr_pions; i++){
			Particle *pion = new Particle(PionType,zero,rcol1);
			if(available_iso-std::abs(iso-iso_system) >= 4){ // pp(pn) pip:0.40(0.3) piz:0.35(0.4) pim:0.25(0.3)
				rdm = Random::shoot();
				if(((iso == 0) && rdm < 0.3) || ((iso == -2) && rdm < 0.40) || rdm < 0.25){
					pion->setType(PiMinus);
					iso_system -= 2;
					available_iso -= 2;
				}
				else if(((iso == 0) && rdm < 0.7) || ((iso == -2) && rdm < 0.75) || rdm < 0.60){
					pion->setType(PiZero);
					available_iso -= 2;
				}
				else{
					pion->setType(PiPlus);
					iso_system += 2;
					available_iso -= 2;
				}
			}
			else if(available_iso-std::abs(iso-iso_system) == 2){
				rdm = Random::shoot();
				//				pn					   		pp too high (nn too low) -> PiMinus(PiPlus)			            			pp too low (nn too high) -> PiPlus(PiMinus)
				if(((iso == 0) && (rdm*0.7 < 0.3)) || ((rdm*0.60 < 0.25) && (Math::sign(iso-iso_system)*2-iso != 0)) || ((rdm*0.75 < 0.40) && (Math::sign(iso-iso_system)*2+iso != 0) && (iso != 0))){
					pion->setType(ParticleTable::getPionType(Math::sign(iso-iso_system)*2));
					iso_system += Math::sign(iso-iso_system)*2;
					available_iso -= 2;
				}
				else{
					PionType = PiZero;
					available_iso -= 2;
				}
			}
			else if(available_iso-std::abs(iso-iso_system) == 0){
				pion->setType(ParticleTable::getPionType(Math::sign(iso-iso_system)*2));
				iso_system += Math::sign(iso-iso_system)*2;
				available_iso -= 2;
			}
			else INCL_ERROR("Pion Generation Problem in NNToMissingStrangenessChannel" << '\n');
			list.push_back(pion);
		}
		
		if(particle2->isLambda()){ // N-K-Lambda
// assert(available_iso == 2);
			if(std::abs(iso-iso_system) == 2){
				particle1->setType(ParticleTable::getNucleonType((iso-iso_system)/2));
				kaon->setType(ParticleTable::getKaonType((iso-iso_system)/2));
			}
			else if(std::abs(iso-iso_system) == 0){
				rdm = G4int(Random::shoot()*2.)*2-1;
				particle1->setType(ParticleTable::getNucleonType(G4int(rdm)));
				kaon->setType(ParticleTable::getKaonType(G4int(-rdm)));
			}
			else INCL_ERROR("Isospin non-conservation in NNToMissingStrangenessChannel" << '\n');
		}
		else if(min_pions == 1){ // N-N-K-Kb chosen
// assert(available_iso == 4);
			Particle *antikaon = new Particle(KMinus,zero,rcol1);
			if(std::abs(iso-iso_system) == 4){
				particle1->setType(ParticleTable::getNucleonType((iso-iso_system)/4));
				particle2->setType(ParticleTable::getNucleonType((iso-iso_system)/4));
				kaon->setType(ParticleTable::getKaonType((iso-iso_system)/4));
				antikaon->setType(ParticleTable::getAntiKaonType((iso-iso_system)/4));
			}
			else if(std::abs(iso-iso_system) == 2){ // choice: kaon type free, nucleon type fixed
				rdm = G4int(Random::shoot()*2.)*2-1;
				particle1->setType(ParticleTable::getNucleonType((iso-iso_system)/2));
				particle2->setType(ParticleTable::getNucleonType((iso-iso_system)/2));
				kaon->setType(ParticleTable::getKaonType(G4int(rdm)));
				antikaon->setType(ParticleTable::getAntiKaonType(G4int(-rdm)));
			}
			else if(std::abs(iso-iso_system) == 0){ // particle1 3/4 proton, 1/4 neutron; particle2 1/4 proton, 3/4 neutron
				rdm = G4int(Random::shoot()*2.)*2-1;
				G4double rdm2 = G4int(Random::shoot()*2.)*2-1;
				particle1->setType(ParticleTable::getNucleonType(std::max(rdm,rdm2)));
				particle2->setType(ParticleTable::getNucleonType(std::min(rdm,rdm2)));
				kaon->setType(ParticleTable::getKaonType(-G4int(rdm)));
				antikaon->setType(ParticleTable::getAntiKaonType(-G4int(rdm2)));
			}
			else INCL_ERROR("Isospin non-conservation in NNToMissingStrangenessChannel" << '\n');
			list.push_back(antikaon);
			nbr_pions += 1; // need for addCreatedParticle loop
		}
		else{// N-K-Sigma
// assert(available_iso == 4);
			if(std::abs(iso-iso_system) == 4){
				particle1->setType(ParticleTable::getNucleonType((iso-iso_system)/4));
				particle2->setType(ParticleTable::getSigmaType((iso-iso_system)/2));
				kaon->setType(ParticleTable::getKaonType((iso-iso_system)/4));
			}
			else if(std::abs(iso-iso_system) == 2){ // choice: sigma quasi-free, kaon type semi-free, nucleon type fixed
				rdm = Random::shoot();
				// pp(pn) sp:0.42(0.23) sz:0.51(0.54) sm:0.07(0.23)
				if(((iso == 0) && (rdm*0.77 < 0.23)) || ((rdm*0.58 < 0.07) && (Math::sign(iso-iso_system)*2-iso != 0)) || ((rdm*0.93 < 0.42) && (Math::sign(iso-iso_system)*2+iso != 0) && (iso != 0))){
					particle2->setType(ParticleTable::getSigmaType(Math::sign(iso-iso_system)*2));
					rdm = G4int(Random::shoot()*2.)*2-1;
					particle1->setType(ParticleTable::getNucleonType(G4int(rdm)));
					kaon->setType(ParticleTable::getKaonType(-G4int(rdm)));
				}
				else{
					particle2->setType(SigmaZero);
					particle1->setType(ParticleTable::getNucleonType((iso-iso_system)/2));
					kaon->setType(ParticleTable::getKaonType((iso-iso_system)/2));
				}
			}
			else if(std::abs(iso-iso_system) == 0){ // choice: sigma free, kaontype semi-free, nucleon type fixed
				if(((iso == 0) && rdm < 0.23) || ((iso == 2) && rdm < 0.42) || rdm < 0.07){
					particle2->setType(SigmaPlus);
					particle1->setType(Neutron);
					kaon->setType(KZero);
				}
				else if(((iso == 0) && rdm < 0.77) || ((iso == 2) && rdm < 0.93) || rdm < 0.58){
					particle2->setType(SigmaZero);
					rdm = G4int(Random::shoot()*2.)*2-1;
					particle1->setType(ParticleTable::getNucleonType(G4int(rdm)));
					kaon->setType(ParticleTable::getKaonType(-G4int(rdm)));
				}
				else{
					particle2->setType(SigmaMinus);
					particle1->setType(Proton);
					kaon->setType(KPlus);
				}
			}
			else INCL_ERROR("Isospin non-conservation in NNToMissingStrangenessChannel" << '\n');
		}
		
		list.push_back(particle1);
		list.push_back(particle2);
		list.push_back(kaon);
		
		PhaseSpaceGenerator::generateBiased(sqrtS, list, list.size()-3, angularSlope);
		
		fs->addModifiedParticle(particle1);
		fs->addModifiedParticle(particle2);
		fs->addCreatedParticle(kaon);
		for(Int_t i=0; i<nbr_pions; i++) fs->addCreatedParticle(list[i]);
		
	}
}
