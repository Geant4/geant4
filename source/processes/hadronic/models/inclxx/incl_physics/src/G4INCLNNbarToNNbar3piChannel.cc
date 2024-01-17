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

#include "G4INCLNNbarToNNbar3piChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	NNbarToNNbar3piChannel::NNbarToNNbar3piChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NNbarToNNbar3piChannel::~NNbarToNNbar3piChannel(){}
	
	void NNbarToNNbar3piChannel::fillFinalState(FinalState *fs) {

		//brief ppbar
        // p pbar -> p pbar pi+ pi- pi0 (BFMM 161)
        // p pbar -> p nbar 2pi- pi+ (BFMM 169)
        // p pbar -> n pbar 2pi+ pi- (BFMM 201)
        // p pbar -> n nbar pi+ pi- pi0 (BFMM 197)
        //
        //brief npbar
        // n pbar -> p pbar 2pi- pi+ (same as BFMM 169)
        // n pbar -> p nbar 2pi- pi0 (same as BFMM 197)
        // n pbar -> n nbar 2pi- pi+ (same as BFMM 169)
        // n pbar -> n pbar pi+ pi- pi0 (same as BFMM 161)
        //
        //brief nnbar
        // n nbar -> n nbar pi+ pi- pi0 (same as BFMM 161)
        // n nbar -> p nbar 2pi- pi+ (same as BFMM 169)
        // n nbar -> n pbar 2pi+ pi- (same as BFMM 201)
        // n nbar -> p pbar pi+ pi- pi0 (same as BFMM 197)
        //
        //brief pnbar
        // p nbar -> p pbar 2pi+ pi- (same as BFMM 169)
        // p nbar -> n pbar 2pi+ pi0 (same as BFMM 197)
        // p nbar -> n nbar 2pi+ pi- (same as BFMM 169)
        // p nbar -> p nbar pi+ pi- pi0 (same as BFMM 161)

		Particle *nucleon;
		Particle *antinucleon;
		
		if(particle1->isNucleon()){
			nucleon = particle1;
			antinucleon = particle2;
		}
		else{
			nucleon = particle2;
			antinucleon = particle1;
		}
		
		const G4double plab = 0.001*KinematicsUtils::momentumInLab(particle1, particle2);
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(nucleon, antinucleon);
		const G4double rdm = Random::shoot();

        const std::vector<G4double> BFMM161 = {-6.434, 1.351, -5.185, 7.754, -1.692, 1.604};
        //const G4double Eth_PPbar_PPbar_pip_pim_pi0 = 1.604;
        const std::vector<G4double> BFMM169 = {3.696, -5.356, -0.053, 1.941, -0.432, 1.624};
        //const G4double Eth_PPbar_PNbar_2pim_pip = 1.624;
        const std::vector<G4double> BFMM201 = {-1.070, -0.636, -0.009, 2.335, -0.499, 1.624};
        //const G4double Eth_PPbar_NPbar_2pip_pim = 1.624;
        const std::vector<G4double> BFMM197 = {1.857, -21.213, -3.448, 0.827, -0.390, 1.616};
        //const G4double Eth_PPbar_NNbar_pip_pim_pi0 = 1.616;

        // pnbar total is same as for npbar
        // ppbar total is same as for nnbar
		const G4double totalppbar = KinematicsUtils::compute_xs(BFMM161, plab) 
		+KinematicsUtils::compute_xs(BFMM169, plab) 
		+KinematicsUtils::compute_xs(BFMM201, plab)
		+KinematicsUtils::compute_xs(BFMM197, plab);
		const G4double totalpnbar = KinematicsUtils::compute_xs(BFMM161, plab) 
		+KinematicsUtils::compute_xs(BFMM197, plab) 
		+2*KinematicsUtils::compute_xs(BFMM169, plab);

		//totalnnbar == totalppbar;
		//totalpnbar == totalnpbar;
		ParticleType Pion1;
		ParticleType Pion2;
		ParticleType Pion3;
		
		//setting types of new particles
		if(nucleon->getType()==Proton){
			if(antinucleon->getType()==antiProton){ // ppbar case
				if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM161, plab)){ // p pbar pi+ pi- pi0 case
					Pion1 = PiMinus;
					Pion2 = PiPlus;
					Pion3 = PiZero;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Proton);
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM161, plab)
				+KinematicsUtils::compute_xs(BFMM169, plab)){ //p nbar 2pi- pi+ case
					Pion1 = PiMinus;
					Pion2 = PiMinus;
					Pion3 = PiPlus;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Proton);
					}
				} 
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM161, plab)
				+KinematicsUtils::compute_xs(BFMM169, plab)
				+KinematicsUtils::compute_xs(BFMM201, plab)){ //n pbar 2pi+ pi- case
					Pion1 = PiPlus;
					Pion2 = PiPlus;
					Pion3 = PiMinus;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Neutron);
					}
				} 
				else{ // n nbar pi+ pi- pi0 case
					Pion1 = PiMinus;
					Pion2 = PiPlus;
					Pion3 = PiZero;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Neutron);
					}
				}
			}
			else{ //antiNeutron (pnbar case)
				if(rdm*totalpnbar < KinematicsUtils::compute_xs(BFMM169, plab)){ // p pbar 2pi+ pi- case
					Pion1 = PiPlus;
					Pion2 = PiPlus;
					Pion3 = PiMinus;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Proton);
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM169, plab)
				+KinematicsUtils::compute_xs(BFMM197, plab)){ // n pbar 2pi+ pi0 case
					Pion1 = PiPlus;
					Pion2 = PiPlus;
					Pion3 = PiZero;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Neutron);
					}
				} 
				else if(rdm*totalppbar < 2*KinematicsUtils::compute_xs(BFMM169, plab)
				+KinematicsUtils::compute_xs(BFMM197, plab)){ // n nbar 2pi+ pi- case
					Pion1 = PiPlus;
					Pion2 = PiPlus;
					Pion3 = PiMinus;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Neutron);
					}
				} 
				else{ // p nbar pi+ pi- pi0 case
					Pion1 = PiMinus;
					Pion2 = PiPlus;
					Pion3 = PiZero;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Proton);
					}
				}
			}
		}
		else{ // neutron
			if(antinucleon->getType()==antiProton){ //npbar case
				if(rdm*totalpnbar < KinematicsUtils::compute_xs(BFMM169, plab)){ // p pbar 2pi- pi+ case
					Pion1 = PiPlus;
					Pion2 = PiMinus;
					Pion3 = PiMinus;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Proton);
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM169, plab)
				+KinematicsUtils::compute_xs(BFMM197, plab)){ // p nbar 2pi- pi0 case
					Pion1 = PiMinus;
					Pion2 = PiMinus;
					Pion3 = PiZero;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Proton);
					}
				} 
				else if(rdm*totalppbar < 2*KinematicsUtils::compute_xs(BFMM169, plab)
				+KinematicsUtils::compute_xs(BFMM197, plab)){ // n nbar 2pi- pi+ case
					Pion1 = PiPlus;
					Pion2 = PiMinus;
					Pion3 = PiMinus;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Neutron);
					}
				} 
				else{ // n pbar pi+ pi- pi0 case
					Pion1 = PiMinus;
					Pion2 = PiPlus;
					Pion3 = PiZero;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Neutron);
					}
				}
			}
			else{ //antiNeutron (nnbar case)
				if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM161, plab)){ // n nbar pi+ pi- pi0 case
					Pion1 = PiMinus;
					Pion2 = PiPlus;
					Pion3 = PiZero;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Neutron);
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM161, plab)
				+KinematicsUtils::compute_xs(BFMM169, plab)){ //p nbar 2pi- pi+ case
					Pion1 = PiMinus;
					Pion2 = PiMinus;
					Pion3 = PiPlus;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Proton);
					}
				} 
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM161, plab)
				+KinematicsUtils::compute_xs(BFMM169, plab)
				+KinematicsUtils::compute_xs(BFMM201, plab)){ //n pbar 2pi+ pi- case
					Pion1 = PiPlus;
					Pion2 = PiPlus;
					Pion3 = PiMinus;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Neutron);
					}
				} 
				else{ // p pbar pi+ pi- pi0 case
					Pion1 = PiMinus;
					Pion2 = PiPlus;
					Pion3 = PiZero;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Proton);
					}
				}
			}
		}
		
		ParticleList list;
		list.push_back(nucleon);
		list.push_back(antinucleon);
		const ThreeVector &rcol = nucleon->getPosition();
		const ThreeVector zero;

		// Create three particle pointers
		Particle *pion1 = nullptr;
		Particle *pion2 = nullptr;
		Particle *pion3 = nullptr;

		// Determine the types of particles based on the random number
		if (rdm < 1.0 / 3.0) {
		    pion1 = new Particle(Pion1, zero, rcol);
		    pion2 = new Particle(Pion2, zero, rcol);
		    pion3 = new Particle(Pion3, zero, rcol);
		} else if (rdm < 2.0 / 3.0) {
		    pion1 = new Particle(Pion1, zero, rcol);
		    pion2 = new Particle(Pion3, zero, rcol);
		    pion3 = new Particle(Pion2, zero, rcol);
		} else {
		    pion1 = new Particle(Pion2, zero, rcol);
		    pion2 = new Particle(Pion1, zero, rcol);
		    pion3 = new Particle(Pion3, zero, rcol);
		}

		list.push_back(pion1);
		list.push_back(pion2);
		list.push_back(pion3);
		
		PhaseSpaceGenerator::generate(sqrtS, list);
		
		fs->addModifiedParticle(nucleon);
		fs->addModifiedParticle(antinucleon);
		fs->addCreatedParticle(pion1);
		fs->addCreatedParticle(pion2);
		fs->addCreatedParticle(pion3);
				
	}
}
