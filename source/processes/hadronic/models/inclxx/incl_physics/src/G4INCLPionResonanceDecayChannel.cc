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

#include "G4INCLPionResonanceDecayChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLPhaseSpaceGenerator.hh"
// #include <cassert>


namespace G4INCL {
	
	const G4double PionResonanceDecayChannel::angularSlope = 0.; // Slope to be defined, if needed.
	
	PionResonanceDecayChannel::PionResonanceDecayChannel(Particle *p, ThreeVector const &dir)
    :theParticle(p), incidentDirection(dir)
	{ }
	
	PionResonanceDecayChannel::~PionResonanceDecayChannel() {}
	
//	Decay during the intranuclear cascade for Omega only
  G4double PionResonanceDecayChannel::computeDecayTime(Particle *p) {
  const G4double m = p->getMass();
  const G4double geff = p->getEnergy()/m;
//		const G4double geta = 1.31e-3;
		const G4double gomega = 8.49;
		G4double gg=0.;
		switch (p->getType()) {
/*			case Eta:
				gg=geta;
				break;*/
			case Omega:
				gg=gomega;
				break;
			default:
				INCL_FATAL("Unrecognized pion resonance type; type=" << p->getType() << '\n');
				break;
		}
		const G4double tpires = -G4INCL::PhysicalConstants::hc/gg*std::log(Random::shoot())*geff;
		return tpires;
	}
	
	void PionResonanceDecayChannel::sampleAngles(G4double *ctet_par, G4double *stet_par, G4double *phi_par) {

		(*ctet_par) = -1.0 + 2.0*Random::shoot();
		if(std::abs(*ctet_par) > 1.0) (*ctet_par) = Math::sign(*ctet_par);
		(*stet_par) = std::sqrt(1.-(*ctet_par)*(*ctet_par));
		(*phi_par) = Math::twoPi * Random::shoot();
	}
	
	void PionResonanceDecayChannel::fillFinalState(FinalState *fs) {
		
		ParticleType createdType;
		ParticleType pionType1=Neutron; // to avoid forgetting pionType definition when 3 particles are emitted
		ParticleType pionType2=Neutron; 

		const G4double sqrtS = theParticle->getMass();
		G4int nbpart = 3; // number of emitted particles
		G4double drnd=Random::shoot();
		switch (theParticle->getType()) {
			case Eta:
				if (drnd < 0.3972) { // renormalized to the only four decays taken into account here
//		2 photons
					nbpart=2;
					theParticle->setType(Photon);
					createdType = Photon;
				}
				else if (drnd < 0.7265) {
//		3 pi0
					theParticle->setType(PiZero);
					pionType1 = PiZero;
					pionType2 = PiZero;
				}
				else if (drnd < 0.9575) {
//		pi+ pi- pi0
				theParticle->setType(PiZero);
				pionType1 = PiPlus;
				pionType2 = PiMinus;
			}
				else {
//		pi+ pi- photon
					theParticle->setType(Photon);
					pionType1 = PiPlus;
					pionType2 = PiMinus;
				}
			break;
			case Omega:
				if (drnd < 0.9009) { // renormalized to the only three decays taken into account here
//		pi+ pi- pi0
					theParticle->setType(PiZero);
					pionType1 = PiPlus;
					pionType2 = PiMinus;
				}
				else if (drnd < 0.9845) {
//		pi0 photon
					nbpart=2;
					theParticle->setType(PiZero);
					createdType = Photon;
				}
				else {
//		pi+ pi-
					nbpart=2;
					theParticle->setType(PiPlus);
					createdType = PiMinus;
				}
				break;
			default:
				INCL_FATAL("Unrecognized pion resonance type; type=" << theParticle->getType() << '\n');
				break;
		}

		if (nbpart == 2) {			
			G4double fi, ctet, stet;
			sampleAngles(&ctet, &stet, &fi);
			
			G4double cfi = std::cos(fi);
			G4double sfi = std::sin(fi);
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
				G4double t1 = ctet+cal*stet*sfi/sal;
				G4double t2 = stet/sal;
				q1=(b1*t1+b2*t2*cfi)/beta;
				q2=(b2*t1-b1*t2*cfi)/beta;
				q3=(b3*t1/beta-t2*sfi);
			} else {
				q1 = stet*cfi;
				q2 = stet*sfi;
				q3 = ctet;
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
		else if (nbpart == 3) {
// assert(pionType1!=Neutron && pionType2!=Neutron);
			ParticleList list;
			list.push_back(theParticle);
			const ThreeVector &rposdecay = theParticle->getPosition();
			const ThreeVector zero;
			Particle *Pion1 = new Particle(pionType1,zero,rposdecay);
			Particle *Pion2 = new Particle(pionType2,zero,rposdecay);
			list.push_back(Pion1);
			list.push_back(Pion2);
			
			fs->addModifiedParticle(theParticle);
			fs->addCreatedParticle(Pion1);
			fs->addCreatedParticle(Pion2);

			//			PhaseSpaceGenerator::generateBiased(sqrtS, list, 0, angularSlope); Biasing?
			PhaseSpaceGenerator::generate(sqrtS, list);
		}

	}
}
