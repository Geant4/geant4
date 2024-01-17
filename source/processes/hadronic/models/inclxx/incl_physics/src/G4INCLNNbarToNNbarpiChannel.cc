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

#include "G4INCLNNbarToNNbarpiChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	NNbarToNNbarpiChannel::NNbarToNNbarpiChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NNbarToNNbarpiChannel::~NNbarToNNbarpiChannel(){}
	
	void NNbarToNNbarpiChannel::fillFinalState(FinalState *fs) {

        //brief ppbar
        // p pbar -> p pbar pi0 (BFMM 185)
        // p pbar -> p nbar pi- (BFMM 188)
        // p pbar -> n pbar pi+ (BFMM 199)
        // p pbar -> n nbar pi0 (no data)
        //
        //brief npbar
        // n pbar -> p pbar pi- (BFMM 491)
        // n pbar -> p nbar pion (impossible)
        // n pbar -> n pbar pi0 (BFMM 495)
        // n pbar -> n nbar pi- (same as BFMM 188)
        //
        //brief nnbar
        // n nbar -> n nbar pi0 (same as BFMM 185)
        // n nbar -> p nbar pi- (same as BFMM 188)
        // n nbar -> n pbar pi+ (same as BFMM 199)
        // n nbar -> p pbar pi0 (no data)
        //
        //brief pnbar 
        // p nbar -> p pbar pi+ (same as BFMM 491)
        // p nbar -> n pbar pion (impossible)
        // p nbar -> p nbar pi0 (same as BFMM 495)
        // p nbar -> n nbar pi+ (same as BFMM 188)

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

        const std::vector<G4double> BFMM185 = {-0.734, 0.841, 0.905, 3.415, -2.316, 0.775};
    //{22.781, -22.602, -0.752, -11.036, 1.548, 0.775};
        //const G4double Eth_PPbar_PPbar_pi0 = 0.775;
        const std::vector<G4double> BFMM188 = { -0.442, 0.501, 0.002, 3.434, -1.201, 0.798};
        //const G4double Eth_PPbar_PNbar_pim = 0.798;
        const std::vector<G4double> BFMM199 = {-2.025, 2.055, -2.355, 6.064, -2.004, 0.798};
        //const G4double Eth_PPbar_NPbar_pip = 0.798;
        const std::vector<G4double> BFMM491 = {24.125, -20.669, -1.534, -19.573, 4.493, 0.787};
        //const G4double Eth_NPbar_PPbar_pim = 0.787; 
        const std::vector<G4double> BFMM495 = {-0.650, -0.140, -0.058, 5.166, -1.705, 0.777};
        //const G4double Eth_NPbar_NPbar_pi0 = 0.777;

        // pnbar total is same as for npbar
        // ppbar total is same as for nnbar
		const G4double totalppbar = KinematicsUtils::compute_xs(BFMM199, plab) +KinematicsUtils::compute_xs(BFMM185, plab) +KinematicsUtils::compute_xs(BFMM188, plab);
		const G4double totalpnbar = KinematicsUtils::compute_xs(BFMM491, plab) +KinematicsUtils::compute_xs(BFMM495, plab) +KinematicsUtils::compute_xs(BFMM188, plab);
		//totalnnbar == totalppbar;
		//totalpnbar == totalnpbar;
		ParticleType PionType;
		
		//setting types of new particles
		if(nucleon->getType()==Proton){
			if(antinucleon->getType()==antiProton){ // ppbar case
				if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM185, plab)){ // ppbarpi0 case
					PionType = PiZero;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Proton);
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM185, plab)+KinematicsUtils::compute_xs(BFMM188, plab)){ //pnbarpi- case
					PionType = PiMinus;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Proton);
					}
				} 
				else{ // npbarpi+ case
					PionType = PiPlus;
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
			else{ //antiNeutron (pnbar case)
				if(rdm*totalpnbar < KinematicsUtils::compute_xs(BFMM491, plab)){ // ppbarpi+ case
					PionType = PiPlus;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Proton);
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM491, plab)+KinematicsUtils::compute_xs(BFMM495, plab)){ //pnbarpi0 case
					PionType = PiZero;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Proton);
					}
				} 
				else{ // nnbarpi+ case
					PionType = PiPlus;
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
		}
		else{ // neutron
			if(antinucleon->getType()==antiProton){ //npbar case
				if(rdm*totalpnbar < KinematicsUtils::compute_xs(BFMM491, plab)){ // ppbarpi- case
					PionType = PiMinus;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Proton);
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM491, plab)+KinematicsUtils::compute_xs(BFMM495, plab)){ //npbarpi0 case
					PionType = PiZero;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiProton);
					}
					else{
						nucleon->setType(antiProton);
						antinucleon->setType(Neutron);
					}
				} 
				else{ // nnbarpi- case
					PionType = PiMinus;
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
			else{ //antiNeutron (nnbar case)
				if(rdm*totalpnbar < KinematicsUtils::compute_xs(BFMM185, plab)){ // nnbarpi0 case
					PionType = PiZero;
					if(rdm<0.5){
						nucleon->setType(Neutron);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Neutron);
					}
				}
				else if(rdm*totalpnbar < KinematicsUtils::compute_xs(BFMM185, plab)+KinematicsUtils::compute_xs(BFMM188, plab)){ //pnbarpi- case
					PionType = PiMinus;
					if(rdm<0.5){
						nucleon->setType(Proton);
						antinucleon->setType(antiNeutron);
					}
					else{
						nucleon->setType(antiNeutron);
						antinucleon->setType(Proton);
					}
				} 
				else{ // npbarpi+ case
					PionType = PiPlus;
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
		}
		
		ParticleList list;
		list.push_back(nucleon);
		list.push_back(antinucleon);
		const ThreeVector &rcol = nucleon->getPosition();
		const ThreeVector zero;
		Particle *pion = new Particle(PionType,zero,rcol);
		list.push_back(pion);
		
		PhaseSpaceGenerator::generate(sqrtS, list);
		
		fs->addModifiedParticle(nucleon);
		fs->addModifiedParticle(antinucleon);
		fs->addCreatedParticle(pion);
				
	}
}
