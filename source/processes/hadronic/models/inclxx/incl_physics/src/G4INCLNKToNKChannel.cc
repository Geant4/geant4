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

#include "G4INCLNKToNKChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	NKToNKChannel::NKToNKChannel(Particle *p1, Particle *p2)
		: particle1(p1), particle2(p2)
		{}
	
	NKToNKChannel::~NKToNKChannel(){}
	
	void NKToNKChannel::fillFinalState(FinalState *fs) {
		
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
		
// assert(ParticleTable::getIsospin(nucleon->getType()) + ParticleTable::getIsospin(kaon->getType()) == 0);
		
		if(kaon->getType() == KZero){
			nucleon->setType(Neutron);
			kaon->setType(KPlus);
		}
		else{
			nucleon->setType(Proton);
			kaon->setType(KZero);
		}
		
		ThreeVector mom_kaon;
		
//		const G4double pLab = KinematicsUtils::momentumInLab(kaon,nucleon);

/*		if(pLab==0) mom_kaon = Random::normVector();
		else{
			const G4double x = kaon->getMomentum().getX();
			const G4double y = kaon->getMomentum().getY();
			const G4double z = kaon->getMomentum().getZ();
			
			const G4double r = std::sqrt(x*x+y*y+z*z);
			const G4double rho = std::sqrt(x*x+y*y);
			
			const G4double b = 12. * pLab/2375.; // correspond to the forward slope description at 2375 MeV/c in K- p elastic
			const G4double cos_theta = std::log(Random::shoot()*(std::exp(b)-std::exp(-b))+std::exp(-b))/b;
			const G4double sin_theta = std::sqrt(1-cos_theta*cos_theta);
			
			const G4double cos_phi = std::cos(Random::shoot()*Math::twoPi);
			const G4double sin_phi = std::sqrt(1-cos_phi*cos_phi);
			
			if(rho == 0) mom_kaon = ThreeVector(sin_theta*cos_phi,sin_theta*sin_phi,cos_theta);
			// Rotation in the direction of the incident kaon
			const G4double px = x/r*cos_theta - y/rho*sin_theta*cos_phi + z/r*x/rho*sin_theta*sin_phi;
			const G4double py = y/r*cos_theta + x/rho*sin_theta*cos_phi + z/r*y/rho*sin_theta*sin_phi;
			const G4double pz = z/r*cos_theta - rho/r*sin_theta*sin_phi;
			
			mom_kaon = ThreeVector(px,py,pz);
		}*/
		
		mom_kaon = Random::normVector();
		
		G4double norm = KinematicsUtils::momentumInCM(kaon,nucleon);

		kaon->setMomentum(mom_kaon*norm);
		nucleon->setMomentum(-mom_kaon*norm);
		
		kaon->adjustEnergyFromMomentum();
		nucleon->adjustEnergyFromMomentum();
		
		fs->addModifiedParticle(nucleon);
		fs->addModifiedParticle(kaon);
				
	}
}
