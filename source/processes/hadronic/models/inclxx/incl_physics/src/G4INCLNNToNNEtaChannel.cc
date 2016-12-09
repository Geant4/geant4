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

#include "G4INCLNNToNNEtaChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	const G4double NNToNNEtaChannel::angularSlope = 6.; // why 6.? same as in multipions...
	
	NNToNNEtaChannel::NNToNNEtaChannel(Particle *p1, Particle *p2)
    : iso1(0),
    iso2(0),
    particle1(p1),
    particle2(p2)
	{}
	
	NNToNNEtaChannel::~NNToNNEtaChannel(){
		
	}
	
	void NNToNNEtaChannel::fillFinalState(FinalState *fs) {
		
		iso1=ParticleTable::getIsospin(particle1->getType());
		iso2=ParticleTable::getIsospin(particle2->getType());
		
		ParticleList list;
		list.push_back(particle1);
		list.push_back(particle2);
		fs->addModifiedParticle(particle1);
		fs->addModifiedParticle(particle2);
				
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(particle1, particle2);

		const ParticleType tn1=ParticleTable::getNucleonType(iso1);
		particle1->setType(tn1);
		const ParticleType tn2=ParticleTable::getNucleonType(iso2);
		particle2->setType(tn2);
		const ThreeVector &rcolnucleon1 = particle1->getPosition();
		const ThreeVector &rcolnucleon2 = particle2->getPosition();
		const ThreeVector rcol = (rcolnucleon1+rcolnucleon2)*0.5;
		const ThreeVector zero;
		ParticleType etaType=Eta;
		Particle *etaCreated = new Particle(etaType,zero,rcol);
		list.push_back(etaCreated);
		fs->addCreatedParticle(etaCreated);

		G4int biasIndex = ((Random::shoot()<0.5) ? 0 : 1);
		PhaseSpaceGenerator::generateBiased(sqrtS, list, biasIndex, angularSlope);
		
	}
		
}
