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

#include "G4INCLNSToNLChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	NSToNLChannel::NSToNLChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NSToNLChannel::~NSToNLChannel(){}
	
	void NSToNLChannel::fillFinalState(FinalState *fs) {
		
		Particle *nucleon;
		Particle *sigma;
		
		if(particle1->isNucleon()){
			nucleon = particle1;
			sigma = particle2;
		}
		else{
			nucleon = particle2;
			sigma = particle1;
		}
		
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(nucleon, sigma);
		
		const G4int iso = ParticleTable::getIsospin(nucleon->getType()) + ParticleTable::getIsospin(sigma->getType());
// assert(iso == -1 || iso == 1);
		
		nucleon->setType(ParticleTable::getNucleonType(iso));
		sigma->setType(Lambda);
		
		G4double mn=nucleon->getMass();
		G4double my=sigma->getMass();
		
		G4double ey=(sqrtS*sqrtS+my*my-mn*mn)/(2*sqrtS);
		G4double en=std::sqrt(ey*ey-my*my+mn*mn);
		nucleon->setEnergy(en);
		sigma->setEnergy(ey);
		G4double py=std::sqrt(ey*ey-my*my);
		
		ThreeVector mom_hyperon = Random::normVector(py);

		sigma->setMomentum(mom_hyperon);
		nucleon->setMomentum(-mom_hyperon);
		
		fs->addModifiedParticle(nucleon);
		fs->addModifiedParticle(sigma);
				
	}
}
