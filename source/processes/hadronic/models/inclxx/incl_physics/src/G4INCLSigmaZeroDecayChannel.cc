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

#include "G4INCLSigmaZeroDecayChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {
	
	SigmaZeroDecayChannel::SigmaZeroDecayChannel(Particle *p, ThreeVector const &dir)
    :theParticle(p), incidentDirection(dir)
	{}
	
	SigmaZeroDecayChannel::~SigmaZeroDecayChannel() {}
	
	
	G4double SigmaZeroDecayChannel::computeDecayTime(Particle *p) {
// assert(p->getType() == SigmaZero);
		const G4double gamma = std::sqrt(1+std::pow(p->getMomentum().mag()/p->getMass(),2));
		const G4double tau = ParticleTable::getWidth(SigmaZero)*3E8*1E15*gamma; // fm
		const G4double t = -tau * std::log(Random::shoot());
		return t;
	}

	void SigmaZeroDecayChannel::sampleAngles(G4double *ctet_par, G4double *stet_par, G4double *phi_par) {

		(*ctet_par) = -1.0 + 2.0*Random::shoot();
		if(std::abs(*ctet_par) > 1.0) (*ctet_par) = Math::sign(*ctet_par); // needed?
		(*stet_par) = std::sqrt(1.-(*ctet_par)*(*ctet_par));
		(*phi_par) = Math::twoPi * Random::shoot();
	}
	
	void SigmaZeroDecayChannel::fillFinalState(FinalState *fs) {
		
// assert( theParticle->getType() == SigmaZero);
		ParticleType createdType = Photon;
		
		const G4double sqrtS = theParticle->getMass();
		
		theParticle->setType(Lambda);
		G4double phi, c_tet, s_tet;
		sampleAngles(&c_tet, &s_tet, &phi);
		
		G4double c_phi = std::cos(phi);
		G4double s_phi = std::sin(phi);
		G4double beta = incidentDirection.mag();
		
		G4double q1, q2, q3;
		G4double sal=0.0;
		if (beta >= 1.0e-10)
			sal = incidentDirection.perp()/beta;
		if (sal >= 1.0e-6) {
			G4double b1 = incidentDirection.getX();
			G4double b2 = incidentDirection.getY();
			G4double b3 = incidentDirection.getZ();
			G4double cal = b3/beta;
			G4double t1 = c_tet+cal*s_tet*s_phi/sal;
			G4double t2 = s_tet/sal;
			q1=(b1*t1+b2*t2*c_phi)/beta;
			q2=(b2*t1-b1*t2*c_phi)/beta;
			q3=(b3*t1/beta-t2*s_phi);
		} else {
			q1 = s_tet*c_phi;
			q2 = s_tet*s_phi;
			q3 = c_tet;
		}
					
		G4double xq = KinematicsUtils::momentumInCM(sqrtS,
												  theParticle->getMass(),
												  ParticleTable::getINCLMass(createdType));
		q1 *= xq;
		q2 *= xq;
		q3 *= xq;
		
		ThreeVector createdMomentum(q1, q2, q3);
		ThreeVector createdPosition(theParticle->getPosition());
		Particle *createdParticle = new Particle(createdType, createdMomentum, createdPosition);
		theParticle->setMomentum(-createdMomentum);
		theParticle->adjustEnergyFromMomentum();
		
		fs->addModifiedParticle(theParticle);
		fs->addCreatedParticle(createdParticle);
		
	}
}
