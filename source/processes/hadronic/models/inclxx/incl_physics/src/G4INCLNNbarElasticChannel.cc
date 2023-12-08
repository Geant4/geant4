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

#include "G4INCLNNbarElasticChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	NNbarElasticChannel::NNbarElasticChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NNbarElasticChannel::~NNbarElasticChannel(){}
	
	void NNbarElasticChannel::fillFinalState(FinalState *fs) {

        //brief ppbar
        // p pbar -> p pbar (BFMM 2)
        //
        //brief npbar
        // n pbar -> n pbar (BFMM 472)
        //
        //brief nnbar
        // n nbar -> n nbar (same as BFMM 2)
        //
        //brief pnbar 
        // p nbar -> p nbar (same as BFMM 472)
        //
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
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(nucleon, antinucleon);
		
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
