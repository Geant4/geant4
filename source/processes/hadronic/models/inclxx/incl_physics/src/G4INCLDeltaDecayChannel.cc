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

#include "G4INCLDeltaDecayChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  DeltaDecayChannel::DeltaDecayChannel(Particle *p, ThreeVector const &dir)
    :theParticle(p), incidentDirection(dir)
  { }

  DeltaDecayChannel::~DeltaDecayChannel() {}

  G4double DeltaDecayChannel::computeDecayTime(Particle *p) {
    const G4double m = p->getMass();
    const G4double g0 = 115.0;
    G4double gg = g0;
    if(m > 1500.0) gg = 200.0;
    const G4double geff = p->getEnergy()/m;
    const G4double qqq = KinematicsUtils::momentumInCM(m, ParticleTable::effectiveNucleonMass, ParticleTable::effectivePionMass);
    const G4double psf = std::pow(qqq, 3)/(std::pow(qqq, 3) + 5832000.0); // phase space factor    5.832E6 = 180^3
    const G4double tdel = -G4INCL::PhysicalConstants::hc/(gg*psf)*std::log(Random::shoot())*geff; // fm
    if( m > 1400) return tdel * 1./(1. + std::pow((m-1400)/g0,2)); // reduction of Delta life time for high masses.
    return tdel; // fm
  }

  void DeltaDecayChannel::sampleAngles(G4double *ctet_par, G4double *stet_par, G4double *phi_par) {
    const G4double hel = theParticle->getHelicity();
    unsigned long loopCounter = 0;
    const unsigned long maxLoopCounter = 10000000;
    do {
      (*ctet_par) = -1.0 + 2.0*Random::shoot();
      if(std::abs(*ctet_par) > 1.0) (*ctet_par) = Math::sign(*ctet_par);
      ++loopCounter;
    } while(loopCounter<maxLoopCounter && Random::shoot() > ((1.0 + 3.0 * hel * (*ctet_par) * (*ctet_par)) /(1.0 + 3.0 * hel))); /* Loop checking, 10.07.2015, D.Mancusi */
    (*stet_par) = std::sqrt(1.-(*ctet_par)*(*ctet_par));
    (*phi_par) = Math::twoPi * Random::shoot();
  }

  void DeltaDecayChannel::fillFinalState(FinalState *fs) {
    //      SUBROUTINE DECAY2(P1,P2,P3,WP,ij,
    //     s       X1,X2,hel,B1,B2,B3)

    // This routine describes the anisotropic decay of a particle of mass
    // xi into 2 particles of masses x1,x2.
    // The anisotropy is supposed to follow a 1+3*hel*(cos(theta))**2
    // law with respect to the direction of the incoming particle.
    // In the input, p1,p2,p3 is the momentum of particle xi.
    // In the output, p1,p2,p3 is the momentum of particle x1 , while
    // q1,q2,q3 is the momentum of particle x2.

    //  COMMON/bl12/QQ1(200),QQ2(200),QQ3(200),QQ4(200),
    // s            YY1(200),YY2(200),YY3(200),YM(200),IPI(200)
    //   common/hazard/ial,IY1,IY2,IY3,IY4,IY5,IY6,IY7,IY8,IY9,IY10,
    // s               IY11,IY12,IY13,IY14,IY15,IY16,IY17,IY18,IY19

    // DATA IY8,IY9,IY10/82345,92345,45681/
    // PCM(E,A,C)=0.5*SQRT((E**2-(A+C)**2)*(E**2-(A-C)**2))/E            P-N20800
    // XI=YM(ij)

    // XE=WP                                                             P-N20810
    // B1=P1/XE                                                          P-N20820
    // B2=P2/XE                                                          P-N20830
    // B3=P3/XE
    // XQ=PCM(XI,X1,X2)

    const G4double deltaMass = theParticle->getMass();

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
    theParticle->setHelicity(0.0);

    ParticleType pionType;
    G4int deltaPDGCode = 0;
    switch(theParticle->getType()) {
      case DeltaPlusPlus:
        theParticle->setType(Proton);
        pionType = PiPlus;
	deltaPDGCode = 2224;
        break;
      case DeltaPlus:
        if(Random::shoot() < 1.0/3.0) {
          theParticle->setType(Neutron);
          pionType = PiPlus;
        } else {
          theParticle->setType(Proton);
          pionType = PiZero;
        }
	deltaPDGCode = 2214;
        break;
      case DeltaZero:
        if(Random::shoot() < 1.0/3.0) {
          theParticle->setType(Proton);
          pionType = PiMinus;
        } else {
          theParticle->setType(Neutron);
          pionType = PiZero;
        }
	deltaPDGCode = 2114;
        break;
      case DeltaMinus:
        theParticle->setType(Neutron);
        pionType = PiMinus;
	deltaPDGCode = 1114;
        break;
      default:
        INCL_FATAL("Unrecognized delta type; type=" << theParticle->getType() << '\n');
        pionType = UnknownParticle;
        break;
    }

    G4double xq = KinematicsUtils::momentumInCM(deltaMass,
        theParticle->getMass(),
        ParticleTable::getINCLMass(pionType));

    q1 *= xq;
    q2 *= xq;
    q3 *= xq;

    ThreeVector pionMomentum(q1, q2, q3);
    ThreeVector pionPosition(theParticle->getPosition());
    Particle *pion = new Particle(pionType, pionMomentum, pionPosition);
    theParticle->setMomentum(-pionMomentum);
    theParticle->adjustEnergyFromMomentum();

    // Set the information about the parent resonance for the two daughters
    // (as unique, integer ID, we take the rounded integer of the resonance mass in keV)
    G4int parentResonanceID = static_cast<G4int>(round(deltaMass/CLHEP::keV));
    pion->setParentResonancePDGCode(deltaPDGCode);
    pion->setParentResonanceID(parentResonanceID);
    theParticle->setParentResonancePDGCode(deltaPDGCode);
    theParticle->setParentResonanceID(parentResonanceID);

    fs->addModifiedParticle(theParticle);
    fs->addCreatedParticle(pion);
    //      call loren(q1,q2,q3,b1,b2,b3,wq)
    //      call loren(p1,p2,p3,b1,b2,b3,wp)
    //      qq1(ij)=q1
    //      qq2(ij)=q2
    //      qq3(ij)=q3
    //      qq4(ij)=wq
    //      ym(ij)=xi
    //      RETURN                                                            P-N21120
    //      END                                                               P-N21130
  }
}
