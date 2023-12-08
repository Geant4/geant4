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

#include "G4INCLNNbarToLLbarChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	NNbarToLLbarChannel::NNbarToLLbarChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NNbarToLLbarChannel::~NNbarToLLbarChannel(){}
	
	void NNbarToLLbarChannel::fillFinalState(FinalState *fs) {
		// this channel include all states with lambdas, sigmas and xis and their antiparticles

        //brief ppbar
        // p pbar -> l lbar (BFMM 121)
        // ppbar -> l lbar pi0 (BFMM 113)
        // ppbar -> splus pim lbar || sminusbar pim l (BFMM 136)
        // ppbar -> sminus pip lbar || splusbar l pip (BFMM 146)
        // ppbar -> sp spbar (BFMM 139)
        // ppbar -> sm smbar (BFMM 149) 
        // ppbar -> szero szerobar (BFMM 144)
        // ppbar -> ximinus ximinusbar (BFMM 101)
        // ppbar -> szero lbar || szerobar l (BFMM 143)
        //
        //
        //brief npbar
        // n pbar -> l lbar pi- (BFMM 487)
        // n pbar -> l sbarplus || lbar sminus (BFMM 488)
        //
        //
        //brief nnbar
        // all same as for ppbar
        //
        //
        //brief pnbar
        // p nbar -> l lbar pi+ (same as BFMM 487)
        // p nbar -> l sbarminus || lbar splus (same as BFMM 488)
        //

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
		// ppbar cross sections

		const std::vector<G4double> BFMM121 = {2.379, -2.738, -1.260, -1.915, 0.430, 1.437};
		//const G4double Eth_PPbar_LLbar = 1.437;
		const std::vector<G4double> BFMM113 = {-0.105, 0.000, -5.099, 0.188, -0.050, 1.820};
		//const G4double Eth_PPbar_LLbar_pi0 = 1.820;
		const std::vector<G4double> BFMM139 = {0.142, -0.291, -1.702, -0.058, 0.001, 1.851};
		//const G4double Eth_PPbar_SpSpbar = 1.851;
		const std::vector<G4double> BFMM149 = {1.855, -2.238, -1.002, -1.279, 0.252, 1.896};
		//const G4double Eth_PPbar_SmSmbar = 1.896;
		const std::vector<G4double> BFMM136 = {1.749, -2.506, -1.222, -1.262, 0.274, 2.042};
		//const G4double Eth_PPbar_SpLbar_pim = 2.042;
		const std::vector<G4double> BFMM146 = {1.037, -1.437, -1.155, -0.709, 0.138, 2.065};
		//const G4double Eth_PPbar_SmLbar_pip = 2.065;
		const std::vector<G4double> BFMM143 = {0.652, -1.006, -1.805, -0.537, 0.121, 1.653};



		//const G4double Eth_PPbar_Szero_Lbar = 1.653;
		//fixed due to limited data
		G4double BFMM144; 
		if(plab > 2.0) BFMM144 = 0.008; //sigmazero sigmazerobar
		else BFMM144 = 0.0;
		G4double BFMM101; 
		if(plab > 2.8) BFMM101 = 0.002; //ximinus ximinusbar
		else BFMM101 = 0.0;

		// npbar cross sections (fixed due to limited data)
		G4double BFMM487;
		if(plab > 2.1) BFMM487 = 0.048; //llbar piminus
		else BFMM487 = 0.0;
		G4double BFMM488;
		if(plab > 2.0) BFMM488 = 0.139; //lsigmaminus +cc
		else BFMM488 = 0.0;

		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(nucleon, antinucleon);
		const G4double totalppbar = KinematicsUtils::compute_xs(BFMM113, plab) 
		+KinematicsUtils::compute_xs(BFMM139, plab) +KinematicsUtils::compute_xs(BFMM136, plab)
		+KinematicsUtils::compute_xs(BFMM146, plab)+KinematicsUtils::compute_xs(BFMM143, plab) 
		+KinematicsUtils::compute_xs(BFMM121, plab)+KinematicsUtils::compute_xs(BFMM149, plab)
		+BFMM144 +BFMM101;
		const G4double totalpnbar = BFMM487 + BFMM488;
		const G4double rdm = Random::shoot();
		
		G4bool thirdparticle = false; //set true if we have pion
		ParticleType PionType;
		//setting types of new particles
		if(nucleon->getType()==Proton){
			if(antinucleon->getType()==antiProton){ //ppbar case
				if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)){ //llbar
					nucleon->setType(Lambda);
					antinucleon->setType(antiLambda);
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144){ //sigmazero sigmazerobar
					nucleon->setType(SigmaZero);
					antinucleon->setType(antiSigmaZero);
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101){ //ximinus ximinusbar
					nucleon->setType(XiMinus);
					antinucleon->setType(antiXiMinus);
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101
				+KinematicsUtils::compute_xs(BFMM113, plab)){ //llbar pi0
					nucleon->setType(Lambda);
					antinucleon->setType(antiLambda);
					thirdparticle = true;
					PionType = PiZero;
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101
				+KinematicsUtils::compute_xs(BFMM113, plab) + KinematicsUtils::compute_xs(BFMM136, plab)){ //splus lbar pim || sminusbar l pim
					G4double rdm2 = Random::shoot();
					if(rdm2 > 0.5){
						nucleon->setType(SigmaPlus);
						antinucleon->setType(antiLambda);
						thirdparticle = true;
						PionType = PiMinus;
					}
					else{
						nucleon->setType(antiSigmaMinus);
						antinucleon->setType(Lambda);
						thirdparticle = true;
						PionType = PiMinus;
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101
				+KinematicsUtils::compute_xs(BFMM113, plab) + KinematicsUtils::compute_xs(BFMM136, plab)
				+KinematicsUtils::compute_xs(BFMM146, plab)){ //sminus lbar pip || splussbar l pip
					G4double rdm2 = Random::shoot();
					if(rdm2 > 0.5){
						nucleon->setType(SigmaMinus);
						antinucleon->setType(antiLambda);
						thirdparticle = true;
						PionType = PiPlus;
					}
					else{
						nucleon->setType(antiSigmaPlus);
						antinucleon->setType(Lambda);
						thirdparticle = true;
						PionType = PiPlus;
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101
				+KinematicsUtils::compute_xs(BFMM113, plab) + KinematicsUtils::compute_xs(BFMM136, plab)
				+KinematicsUtils::compute_xs(BFMM146, plab)+KinematicsUtils::compute_xs(BFMM143, plab)){ //szero lbar || szerobar l 
					G4double rdm2 = Random::shoot();
					if(rdm2 > 0.5){
						nucleon->setType(SigmaZero);
						antinucleon->setType(antiLambda);
					}
					else{
						nucleon->setType(antiSigmaZero);
						antinucleon->setType(Lambda);
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101
				+KinematicsUtils::compute_xs(BFMM113, plab) + KinematicsUtils::compute_xs(BFMM136, plab)
				+KinematicsUtils::compute_xs(BFMM146, plab)+KinematicsUtils::compute_xs(BFMM143, plab)
				+KinematicsUtils::compute_xs(BFMM139, plab)){ //sp spbar 
					nucleon->setType(SigmaPlus);
					antinucleon->setType(antiSigmaPlus);
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101
				+KinematicsUtils::compute_xs(BFMM113, plab) + KinematicsUtils::compute_xs(BFMM136, plab)
				+KinematicsUtils::compute_xs(BFMM146, plab)+KinematicsUtils::compute_xs(BFMM143, plab)
				+KinematicsUtils::compute_xs(BFMM139, plab)+KinematicsUtils::compute_xs(BFMM149, plab)){ //sm smbar
					nucleon->setType(SigmaMinus);
					antinucleon->setType(antiSigmaMinus);
				}
				else{
					INCL_ERROR("out of total ppbar sum in LLbar channel");
				}
			}
			else{ //pnbar case charge +1
				if(rdm*totalpnbar < BFMM488){
					G4double rdm2 = Random::shoot();
					if(rdm2 > 0.5){
						nucleon->setType(Lambda);
						antinucleon->setType(antiSigmaMinus); //charge +1
					}
					else{
						nucleon->setType(antiLambda);
						antinucleon->setType(SigmaPlus); //charge +1
					}
				}
				else{
					nucleon->setType(Lambda);
					antinucleon->setType(antiLambda);
					thirdparticle = true;
					PionType = PiPlus;
				}
			}
		}
		else{ // neutron
			if(antinucleon->getType()==antiNeutron){ //nnbar case same as ppbar
				if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)){ //llbar
					nucleon->setType(Lambda);
					antinucleon->setType(antiLambda);
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144){ //sigmazero sigmazerobar
					nucleon->setType(SigmaZero);
					antinucleon->setType(antiSigmaZero);
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101){ //ximinus ximinusbar
					nucleon->setType(XiMinus);
					antinucleon->setType(antiXiMinus);
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101
				+KinematicsUtils::compute_xs(BFMM113, plab)){ //llbar pi0
					nucleon->setType(Lambda);
					antinucleon->setType(antiLambda);
					thirdparticle = true;
					PionType = PiZero;
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101
				+KinematicsUtils::compute_xs(BFMM113, plab) + KinematicsUtils::compute_xs(BFMM136, plab)){ //splus lbar pim || sminusbar l pim
					G4double rdm2 = Random::shoot();
					if(rdm2 > 0.5){
						nucleon->setType(SigmaPlus);
						antinucleon->setType(antiLambda);
						thirdparticle = true;
						PionType = PiMinus;
					}
					else{
						nucleon->setType(antiSigmaMinus);
						antinucleon->setType(Lambda);
						thirdparticle = true;
						PionType = PiMinus;
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101
				+KinematicsUtils::compute_xs(BFMM113, plab) + KinematicsUtils::compute_xs(BFMM136, plab)
				+KinematicsUtils::compute_xs(BFMM146, plab)){ //sminus lbar pip || splussbar l pip
					G4double rdm2 = Random::shoot();
					if(rdm2 > 0.5){
						nucleon->setType(SigmaMinus); //charge -1
						antinucleon->setType(antiLambda);
						thirdparticle = true;
						PionType = PiPlus;
					}
					else{
						nucleon->setType(antiSigmaPlus); //charge -1
						antinucleon->setType(Lambda);
						thirdparticle = true;
						PionType = PiPlus;
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101
				+KinematicsUtils::compute_xs(BFMM113, plab) + KinematicsUtils::compute_xs(BFMM136, plab)
				+KinematicsUtils::compute_xs(BFMM146, plab)+KinematicsUtils::compute_xs(BFMM143, plab)){ //szero lbar || szerobar l 
					G4double rdm2 = Random::shoot();
					if(rdm2 > 0.5){
						nucleon->setType(SigmaZero);
						antinucleon->setType(antiLambda);
					}
					else{
						nucleon->setType(antiSigmaZero);
						antinucleon->setType(Lambda);
					}
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101
				+KinematicsUtils::compute_xs(BFMM113, plab) + KinematicsUtils::compute_xs(BFMM136, plab)
				+KinematicsUtils::compute_xs(BFMM146, plab)+KinematicsUtils::compute_xs(BFMM143, plab)
				+KinematicsUtils::compute_xs(BFMM139, plab)){ //sp spbar 
					nucleon->setType(SigmaPlus);
					antinucleon->setType(antiSigmaPlus);
				}
				else if(rdm*totalppbar < KinematicsUtils::compute_xs(BFMM121, plab)+BFMM144+BFMM101
				+KinematicsUtils::compute_xs(BFMM113, plab) + KinematicsUtils::compute_xs(BFMM136, plab)
				+KinematicsUtils::compute_xs(BFMM146, plab)+KinematicsUtils::compute_xs(BFMM143, plab)
				+KinematicsUtils::compute_xs(BFMM139, plab)+KinematicsUtils::compute_xs(BFMM149, plab)){ //sm smbar
					nucleon->setType(SigmaMinus);
					antinucleon->setType(antiSigmaMinus);
				}
				else{
					INCL_ERROR("out of total nnbar sum in LLbar channel");
				}
			}
			else{ //npbar case charge -1
				if(rdm*totalpnbar < BFMM488){
					G4double rdm2 = Random::shoot();
					if(rdm2 > 0.5){
						nucleon->setType(Lambda);
						antinucleon->setType(antiSigmaPlus); //charge -1
					}
					else{
						nucleon->setType(antiLambda);
						antinucleon->setType(SigmaMinus); //charge -1
					}
				}
				else{
					nucleon->setType(Lambda);
					antinucleon->setType(antiLambda);
					thirdparticle = true;
					PionType = PiMinus;
				}
			}
		}

		//now assigning momentum to the final particles

		if(thirdparticle){ //three particles
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
		else{//only two particles
			G4double mn=nucleon->getMass();
			G4double my=antinucleon->getMass();
			
			G4double ey=(sqrtS*sqrtS+my*my-mn*mn)/(2*sqrtS);
			G4double en=std::sqrt(ey*ey-my*my+mn*mn);
			nucleon->setEnergy(en);
			antinucleon->setEnergy(ey);
			G4double py=std::sqrt(ey*ey-my*my);
			
			ThreeVector mom_antinucleon = Random::normVector(py);

			antinucleon->setMomentum(mom_antinucleon);
			nucleon->setMomentum(-mom_antinucleon);

			fs->addModifiedParticle(nucleon);
			fs->addModifiedParticle(antinucleon);
		}
				
	}
}
