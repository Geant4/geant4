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

#include "G4INCLDeltaProductionChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"

namespace G4INCL {

  const G4int DeltaProductionChannel::maxTries = 100000;

  DeltaProductionChannel::DeltaProductionChannel(Particle *p1,
						 Particle *p2)
    : particle1(p1), particle2(p2)
  {}

  DeltaProductionChannel::~DeltaProductionChannel() {}

  G4double DeltaProductionChannel::sampleDeltaMass(G4double ecm) {
    const G4double maxDeltaMass = ecm - ParticleTable::effectiveNucleonMass - 1.0;
    const G4double maxDeltaMassRndm = std::atan((maxDeltaMass-ParticleTable::effectiveDeltaMass)*2./ParticleTable::effectiveDeltaWidth);
    const G4double deltaMassRndmRange = maxDeltaMassRndm - ParticleTable::minDeltaMassRndm;
// assert(deltaMassRndmRange>0.);

    G4double y=ecm*ecm;
    G4double q2=(y-1.157776E6)*(y-6.4E5)/y/4.0; // 1.157776E6 = 1076^2, 6.4E5 = 800^2
    G4double q3=std::pow(std::sqrt(q2), 3.);
    const G4double f3max=q3/(q3+5.832E6); // 5.832E6 = 180^3
    G4double x;

    G4int nTries = 0;
    G4bool success = false;
    while(!success) { /* Loop checking, 10.07.2015, D.Mancusi */
      if(++nTries >= maxTries) {
        INCL_WARN("DeltaProductionChannel::sampleDeltaMass loop was stopped because maximum number of tries was reached. Minimum delta mass "
                  << ParticleTable::minDeltaMass << " MeV with CM energy " << ecm << " MeV may be unphysical." << '\n');
        return ParticleTable::minDeltaMass;
      }

      G4double rndm = ParticleTable::minDeltaMassRndm + Random::shoot() * deltaMassRndmRange;
      y = std::tan(rndm);
      x = ParticleTable::effectiveDeltaMass + 0.5*ParticleTable::effectiveDeltaWidth*y;
// assert(x>=ParticleTable::minDeltaMass && ecm >= x + ParticleTable::effectiveNucleonMass + 1.0);

      // generation of the delta mass with the penetration factor
      // (see prc56(1997)2431)
      y=x*x;
      q2=(y-1.157776E6)*(y-6.4E5)/y/4.0; // 1.157776E6 = 1076^2, 6.4E5 = 800^2
      q3=std::pow(std::sqrt(q2), 3.);
      const G4double f3=q3/(q3+5.832E6); // 5.832E6 = 180^3
      rndm = Random::shoot();
      if (rndm*f3max < f3)
        success = true;
    }
    return x;
  }

  void DeltaProductionChannel::fillFinalState(FinalState *fs) {
    /**
     * Delta production
     *
     * The production is not isotropic in this version it has the same
     * exp(b*t) structure as the nn elastic scattering (formula 2.3 of
     * j.cugnon et al, nucl phys a352(1981)505) parametrization of b
     * taken from ref. prc56(1997)2431
    */
    //    100 IF (K4.NE.1) GO TO 101 // ThA K4 = 2 by default
    //    ParticleType p1TypeOld = particle1->getType();
    //    ParticleType p2TypeOld = particle2->getType();
    G4double ecm = KinematicsUtils::totalEnergyInCM(particle1, particle2);

    const G4int isospin = ParticleTable::getIsospin(particle1->getType()) +
      ParticleTable::getIsospin(particle2->getType());

    // Calculate the outcome of the channel:
    const ThreeVector &particle1Momentum = particle1->getMomentum();
    G4double pin = particle1Momentum.mag();
    G4double rndm = 0.0, b = 0.0;

    G4double xmdel = sampleDeltaMass(ecm);
    //  deltaProduction103: // This label is not used
    G4double pnorm = KinematicsUtils::momentumInCM(ecm, ParticleTable::effectiveNucleonMass, xmdel);
    if (pnorm <= 0.0) pnorm=0.000001;
    G4int index=0;
    G4int index2=0;
    rndm = Random::shoot();
    if (rndm < 0.5) index=1;
    if (isospin == 0) { // pn case
      rndm = Random::shoot();
      if (rndm < 0.5) index2=1;
    }

    //    G4double x=0.001*0.5*ecm*std::sqrt(ecm*ecm-4.*ParticleTable::effectiveNucleonMass2)
    //      / ParticleTable::effectiveNucleonMass;
    G4double x = 0.001 * KinematicsUtils::momentumInLab(ecm*ecm, ParticleTable::effectiveNucleonMass, ParticleTable::effectiveNucleonMass);
    if(x < 1.4) {
      b=(5.287/(1.+std::exp((1.3-x)/0.05)))*1.e-6;
    } else {
      b=(4.65+0.706*(x-1.4))*1.e-6;
    }
    G4double xkh = 2.*b*pin*pnorm;
    rndm = Random::shoot();
    G4double ctet=1.0+std::log(1.-rndm*(1.-std::exp(-2.*xkh)))/xkh;
    if(std::abs(ctet) > 1.0) ctet = Math::sign(ctet);
    G4double stet = std::sqrt(1.-ctet*ctet);

    rndm = Random::shoot();
    G4double fi = Math::twoPi*rndm;
    G4double cfi = std::cos(fi);
    G4double sfi = std::sin(fi);
    // delta production: correction of the angular distribution 02/09/02

    G4double xx = particle1Momentum.perp2();
    const G4double particle1MomentumZ = particle1Momentum.getZ();
    G4double zz = std::pow(particle1MomentumZ, 2);
    G4double xp1, xp2, xp3;
    if (xx >= zz*1.e-8) {
      G4double yn = std::sqrt(xx);
      G4double zn = yn*pin;
      G4double ex[3], ey[3], ez[3];
      G4double p1 = particle1Momentum.getX();
      G4double p2 = particle1Momentum.getY();
      G4double p3 = particle1MomentumZ;
      ez[0] = p1/pin;
      ez[1] = p2/pin;
      ez[2] = p3/pin;
      ex[0] = p2/yn;
      ex[1] = -p1/yn;
      ex[2] = 0.0;
      ey[0] = p1*p3/zn;
      ey[1] = p2*p3/zn;
      ey[2] = -xx/zn;
      xp1 = (ex[0]*cfi*stet+ey[0]*sfi*stet+ez[0]*ctet)*pnorm;
      xp2 = (ex[1]*cfi*stet+ey[1]*sfi*stet+ez[1]*ctet)*pnorm;
      xp3 = (ex[2]*cfi*stet+ey[2]*sfi*stet+ez[2]*ctet)*pnorm;
    }else {
      xp1=pnorm*stet*cfi;
      xp2=pnorm*stet*sfi;
      xp3=pnorm*ctet;
    }
    // end of correction angular distribution of delta production
    G4double e3 = std::sqrt(xp1*xp1+xp2*xp2+xp3*xp3
			  +ParticleTable::effectiveNucleonMass2);
    //      if(k4.ne.0) go to 161

    // long-lived delta
    if (index != 1) {
      ThreeVector mom(xp1, xp2, xp3);
      particle1->setMomentum(mom);
      //       e1=ecm-eout1
    } else {
      ThreeVector mom(-xp1, -xp2, -xp3);
      particle1->setMomentum(mom);
      //      e1=ecm-eout1
    }

    particle1->setEnergy(ecm - e3);
    particle2->setEnergy(e3);
    particle2->setMomentum(-particle1->getMomentum());

    // SYMMETRIZATION OF CHARGES IN pn -> N DELTA
    // THE TEST ON "INDEX" ABOVE SYMETRIZES THE EXCITATION OF ONE
    // OF THE NUCLEONS WITH RESPECT TO THE DELTA EXCITATION
    // (SEE NOTE 16/10/97)
    G4int is1 = ParticleTable::getIsospin(particle1->getType());
    G4int is2 = ParticleTable::getIsospin(particle2->getType());
    if (isospin == 0) {
      if(index2 == 1) {
        G4int isi=is1;
        is1=is2;
        is2=isi;
      }
      particle1->setHelicity(0.0);
    } else {
      rndm = Random::shoot();
      if (rndm >= 0.25) {
        is1=3*is1;
        is2=-is2;
      }
      particle1->setHelicity(ctet*ctet);
    }

    if(is1 == ParticleTable::getIsospin(DeltaMinus)) {
      particle1->setType(DeltaMinus);
    } else if(is1 == ParticleTable::getIsospin(DeltaZero)) {
      particle1->setType(DeltaZero);
    } else if(is1 == ParticleTable::getIsospin(DeltaPlus)) {
      particle1->setType(DeltaPlus);
    } else if(is1 == ParticleTable::getIsospin(DeltaPlusPlus)) {
      particle1->setType(DeltaPlusPlus);
    }

    if(is2 == ParticleTable::getIsospin(Proton)) {
      particle2->setType(Proton);
    } else if(is2 == ParticleTable::getIsospin(Neutron)) {
      particle2->setType(Neutron);
    }

    if(particle1->isDelta()) particle1->setMass(xmdel);
    if(particle2->isDelta()) particle2->setMass(xmdel);

    fs->addModifiedParticle(particle1);
    fs->addModifiedParticle(particle2);
  }
}
