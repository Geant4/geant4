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

#include "G4INCLNNbarToNNbar2piChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	NNbarToNNbar2piChannel::NNbarToNNbar2piChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NNbarToNNbar2piChannel::~NNbarToNNbar2piChannel(){}
	
	void NNbarToNNbar2piChannel::fillFinalState(FinalState *fs) {

        //brief ppbar
        // p pbar -> p pbar pi+ pi- (BFMM 167)
        // p pbar -> p nbar pi- pi0 (same as BFMM 490)
        // p pbar -> n pbar pi+ pi0 (same as BFMM 490)
        // p pbar -> n nbar pi+ pi- (BFMM 198)
        //
        //brief npbar
        // n pbar -> p pbar pi- pi0 (BFMM 490)
        // n pbar -> p nbar pi- pi- (BFMM 492)
        // n pbar -> n nbar pi- pi0 (same as BFMM 490)
        // n pbar -> n pbar pi+ pi- (BFMM 494)
        //
        //brief nnbar
        // n nbar -> n nbar pi+ pi- (same as BFMM 167)
        // n nbar -> p nbar pi- pi0 (same as BFMM 490)
        // n nbar -> n pbar pi+ pi0 (same as BFMM 490)
        // n nbar -> p pbar pi+ pi- (same as BFMM 198)
        //
        //brief pnbar
        // p nbar -> p pbar pi+ pi0 (same as BFMM 490)
        // p nbar -> n pbar pi+ pi+ (same as BFMM 492)
        // p nbar -> n nbar pi+ pi0 (same as BFMM 490)
        // p nbar -> p nbar pi+ pi- (same as BFMM 494)

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
		
		const G4double plab = 0.001*KinematicsUtils::momentumInLab(particle1, particle2); //GeV
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(nucleon, antinucleon);
		const G4double rdm = Random::shoot();

        const std::vector<G4double> BFMM167 = {-6.885, 0.476, 1.206, 13.857, -5.728, 1.220};
        //const G4double Eth_PPbar_PPbar_pip_pim = 1.220;
        const std::vector<G4double> BFMM198 = {1.857, -21.213, -3.448, 0.827, -0.390, 1.231};
        //const G4double Eth_PPbar_NNbar_pip_pim = 1.231;
        const std::vector<G4double> BFMM490 = {-3.594, 0.811, 0.306, 5.108, -1.625, 1.201};
        //const G4double Eth_PNbar_PPbar_pim_pi0 = 1.201;
        const std::vector<G4double> BFMM492 = {-5.443, 7.254, -2.936, 8.441, -2.588, 1.221};
        //const G4double Eth_PNbar_NPbar_pim_pim = 1.221;
        const std::vector<G4double> BFMM494 = {21.688, -38.709, -2.062, -17.783, 3.895, 1.221};
        //const G4double Eth_NPbar_NPbar_pip_pim = 1.221; 

        // pnbar total is same as for npbar
        // ppbar total is same as for nnbar
		const G4double totalppbar = KinematicsUtils::compute_xs(BFMM167, plab) +KinematicsUtils::compute_xs(BFMM198, plab) +2*KinematicsUtils::compute_xs(BFMM490, plab);
		const G4double totalpnbar = KinematicsUtils::compute_xs(BFMM492, plab) +KinematicsUtils::compute_xs(BFMM494, plab) +2*KinematicsUtils::compute_xs(BFMM490, plab);
		//totalnnbar == totalppbar;
		//totalpnbar == totalnpbar;
		ParticleType Pion1;
		ParticleType Pion2;
		
		//setting types of new particles
		if(nucleon->getType()==Proton){
			if(antinucleon->getType()==antiProton){ // ppbar case
				if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM167, plab)){ // ppbarpi-pi+ case
					Pion1 = PiMinus;
					Pion2 = PiPlus;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Proton);
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM167, plab)+KinematicsUtils::compute_xs(BFMM490, plab)){ //pnbarpi-pi0 case
					Pion1 = PiMinus;
					Pion2 = PiZero;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Proton);
					}
				} 
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM167, plab)+2*KinematicsUtils::compute_xs(BFMM490, plab)){ //npbarpi+pi0 case
					Pion1 = PiPlus;
					Pion2 = PiZero;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Neutron);
					}
				} 
				else{ // n nbar pi+ pi- case case
					Pion1 = PiMinus;
					Pion2 = PiPlus;
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
				if(rdm*totalpnbar < KinematicsUtils::compute_xs(BFMM490, plab)){ // p pbar pi+ pi0 case
					Pion1 = PiZero;
					Pion2 = PiPlus;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Proton);
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM490, plab)+KinematicsUtils::compute_xs(BFMM492, plab)){ // n pbar pi+ pi+ case
					Pion1 = PiPlus;
					Pion2 = PiPlus;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Neutron);
					}
				} 
				else if(rdm*totalppbar < 2*KinematicsUtils::compute_xs(BFMM490, plab)+KinematicsUtils::compute_xs(BFMM492, plab)){ // n nbar pi+ pi0 case
					Pion1 = PiZero;
					Pion2 = PiPlus;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Neutron);
					}
				} 
				else{ // p nbar pi+ pi- case
					Pion1 = PiMinus;
					Pion2 = PiPlus;
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
				if(rdm*totalpnbar < KinematicsUtils::compute_xs(BFMM490, plab)){ // p pbar pi- pi0 case
					Pion1 = PiZero;
					Pion2 = PiMinus;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Proton);
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM490, plab)+KinematicsUtils::compute_xs(BFMM492, plab)){ // p nbar pi- pi- case
					Pion1 = PiMinus;
					Pion2 = PiMinus;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Proton);
					}
				} 
				else if(rdm*totalppbar < 2*KinematicsUtils::compute_xs(BFMM490, plab)+KinematicsUtils::compute_xs(BFMM492, plab)){ // n nbar pi- pi0 case
					Pion1 = PiZero;
					Pion2 = PiMinus;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Neutron);
					}
				} 
				else{ // n pbar pi+ pi- case
					Pion1 = PiMinus;
					Pion2 = PiPlus;
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
				if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM167, plab)){ // nnbarpi-pi+ case
					Pion1 = PiMinus;
					Pion2 = PiPlus;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Neutron);
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM167, plab)+KinematicsUtils::compute_xs(BFMM490, plab)){ //pnbarpi-pi0 case
					Pion1 = PiMinus;
					Pion2 = PiZero;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Proton);
					}
				} 
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM167, plab)+2*KinematicsUtils::compute_xs(BFMM490, plab)){ //npbarpi+pi0 case
					Pion1 = PiPlus;
					Pion2 = PiZero;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Neutron);
					}
				} 
				else{ // p pbar pi+ pi- case
					Pion1 = PiMinus;
					Pion2 = PiPlus;
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

		Particle *pion2 = new Particle(Pion1,zero,rcol);
		Particle *pion1 = new Particle(Pion2,zero,rcol);
		if(rdm < 0.5){
			pion1->setType(Pion1);
			pion2->setType(Pion2);
		}

		list.push_back(pion1);
		list.push_back(pion2);
		
		PhaseSpaceGenerator::generate(sqrtS, list);
		
		fs->addModifiedParticle(nucleon);
		fs->addModifiedParticle(antinucleon);
		fs->addCreatedParticle(pion1);
		fs->addCreatedParticle(pion2);
				
	}
}
