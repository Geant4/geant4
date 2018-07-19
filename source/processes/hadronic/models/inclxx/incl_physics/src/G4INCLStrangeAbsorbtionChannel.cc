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

#include "G4INCLStrangeAbsorbtionChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"

namespace G4INCL {

  StrangeAbsorbtionChannel::StrangeAbsorbtionChannel(Particle *p1, Particle *p2)
    : particle1(p1), particle2(p2)
  {

  }

  StrangeAbsorbtionChannel::~StrangeAbsorbtionChannel(){}
  
  void StrangeAbsorbtionChannel::sampleAngles(G4double *ctet_par, G4double *stet_par, G4double *phi_par) {
	  (*ctet_par) = -1.0 + 2.0*Random::shoot();
	  if(std::abs(*ctet_par) > 1.0) (*ctet_par) = Math::sign(*ctet_par); // needed?
	  (*stet_par) = std::sqrt(1.-(*ctet_par)*(*ctet_par));
	  (*phi_par) = Math::twoPi * Random::shoot();
  }

  void StrangeAbsorbtionChannel::fillFinalState(FinalState *fs) {
    Particle * nucleon;
    Particle * strange;
    
    ThreeVector const incidentDirection = particle1->getMomentum() + particle2->getMomentum();
    
    if(particle1->isNucleon()) {
      nucleon = particle1;
      strange = particle2;
    } else {
      nucleon = particle2;
      strange = particle1;
    }
// assert(strange->isSigma() || strange->isAntiKaon());

    ParticleType finalType = Neutron;
    if(ParticleConfig::isPair(nucleon, strange, Neutron, KZeroBar)) {
      finalType = PiZero;
    } else if(ParticleConfig::isPair(nucleon, strange, Proton, KZeroBar)) {
      finalType = PiPlus;
    } else if(ParticleConfig::isPair(nucleon, strange, Neutron, KMinus)) {
      finalType = PiMinus;
    } else if(ParticleConfig::isPair(nucleon, strange, Proton, KMinus)) {
      finalType = PiZero;
    } else if(ParticleConfig::isPair(nucleon, strange, Proton, SigmaMinus)) {
      finalType = Neutron;
    } else if(ParticleConfig::isPair(nucleon, strange, Neutron, SigmaZero)) {
      finalType = Neutron;
    } else if(ParticleConfig::isPair(nucleon, strange, Proton, SigmaZero)) {
      finalType = Proton;
    } else if(ParticleConfig::isPair(nucleon, strange, Neutron, SigmaPlus)) {
      finalType = Proton;
    } else {
      INCL_ERROR("Unknown particle pair in Strange-N absorbtion: " << nucleon << '\t' << strange << '\n');
      return;
    }

    G4double energycm = KinematicsUtils::totalEnergyInCM(nucleon, strange);
    
    G4double finalTypemass = ParticleTable::getINCLMass(finalType);
    nucleon->setType(Lambda); // nucleon becomes the lambda
    G4double lambdamass = nucleon->getMass();
    
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
	
	G4double xq = KinematicsUtils::momentumInCM(energycm,
											  lambdamass,
											  finalTypemass);
	
	q1 *= xq;
	q2 *= xq;
	q3 *= xq;
	
	ThreeVector finalMomentum(q1, q2, q3);
	
	strange->setType(finalType);
	strange->setMomentum(finalMomentum);
	strange->adjustEnergyFromMomentum();
	nucleon->setMomentum(-finalMomentum);
	nucleon->adjustEnergyFromMomentum();

    fs->addModifiedParticle(nucleon); // nucleon became a lambda
    fs->addModifiedParticle(strange);  // the strange particle became an unstange particle
  }

}
