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
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"

namespace G4INCL {

  void KinematicsUtils::transformToLocalEnergyFrame(Nucleus const * const n, Particle * const p) {
    const G4double localEnergy = KinematicsUtils::getLocalEnergy(n, p);
    const G4double localTotalEnergy = p->getEnergy() - localEnergy;
    p->setEnergy(localTotalEnergy);
    p->adjustMomentumFromEnergy();
  }

  G4double KinematicsUtils::getLocalEnergy(Nucleus const * const n, Particle * const p) {
// assert(!p->isPion()); // No local energy for pions

    G4double vloc = 0.0;
    const G4double r = p->getPosition().mag();
    const G4double mass = p->getMass();

    // Local energy is constant outside the surface
    if(r > n->getUniverseRadius()) {
      WARN("Tried to evaluate local energy for a particle outside the maximum radius."
            << std::endl << p->print() << std::endl
            << "Maximum radius = " << n->getDensity()->getMaximumRadius() << std::endl
            << "Universe radius = " << n->getUniverseRadius() << std::endl);
      return 0.0;
    }

    G4double pfl0 = 0.0;
    const G4double kinE = p->getKineticEnergy();
    if(kinE <= n->getPotential()->getFermiEnergy(p->getType())) {
      pfl0 = n->getPotential()->getFermiMomentum(p);
    } else {
      const G4double tf0 = p->getPotentialEnergy() - n->getPotential()->getSeparationEnergy(p);
      if(tf0<0.0) return 0.0;
      pfl0 = std::sqrt(tf0*(tf0 + 2.0*mass));
    }
    const G4double pl = pfl0*n->getDensity()->getMaxTFromR(r);
    vloc = std::sqrt(pl*pl + mass*mass) - mass;

    return vloc;
  }

  ThreeVector KinematicsUtils::makeBoostVector(Particle const * const p1, Particle const * const p2){
    const G4double totalEnergy = p1->getEnergy() + p2->getEnergy();
    return ((p1->getMomentum() + p2->getMomentum())/totalEnergy);
  }

  G4double KinematicsUtils::totalEnergyInCM(Particle const * const p1, Particle const * const p2){
    return std::sqrt(squareTotalEnergyInCM(p1,p2));
  }

  G4double KinematicsUtils::squareTotalEnergyInCM(Particle const * const p1, Particle const * const p2) {
    G4double beta2 = KinematicsUtils::makeBoostVector(p1, p2).mag2();
    if(beta2 > 1.0) {
      ERROR("KinematicsUtils::squareTotalEnergyInCM: beta2 == " << beta2 << " > 1.0" << std::endl);
      beta2 = 0.0;
    }
    return (1.0 - beta2)*std::pow(p1->getEnergy() + p2->getEnergy(), 2);
  }

  G4double KinematicsUtils::momentumInCM(Particle const * const p1, Particle const * const p2) {
    const G4double m1sq = std::pow(p1->getMass(),2);
    const G4double m2sq = std::pow(p2->getMass(),2);
    const G4double z = p1->getEnergy()*p2->getEnergy() - p1->getMomentum().dot(p2->getMomentum());
    G4double pcm2 = (z*z-m1sq*m2sq)/(2*z+m1sq+m2sq);
    if(pcm2 < 0.0) {
      ERROR("KinematicsUtils::momentumInCM: pcm2 == " << pcm2 << " < 0.0" << std::endl);
      pcm2 = 0.0;
    }
    return std::sqrt(pcm2);
  }

  G4double KinematicsUtils::momentumInCM(const G4double E, const G4double M1, const G4double M2) {
    return 0.5*std::sqrt((E*E - std::pow(M1 + M2, 2))
			 *(E*E - std::pow(M1 - M2, 2)))/E;
  }

  G4double KinematicsUtils::momentumInLab(const G4double s, const G4double m1, const G4double m2) {
    const G4double m1sq = m1*m1;
    const G4double m2sq = m2*m2;
    G4double plab2 = (s*s-2*s*(m1sq+m2sq)+(m1sq-m2sq)*(m1sq-m2sq))/(4*m2sq);
    if(plab2 < 0.0) {
      ERROR("KinematicsUtils::momentumInLab: plab2 == " << plab2 << " < 0.0; m1sq == " << m1sq << "; m2sq == " << m2sq << "; s == " << s << std::endl);
      plab2 = 0.0;
    }
    return std::sqrt(plab2);
  }

  G4double KinematicsUtils::momentumInLab(Particle const * const p1, Particle const * const p2) {
    const G4double m1 = p1->getMass();
    const G4double m2 = p2->getMass();
    const G4double s = squareTotalEnergyInCM(p1, p2);
    return KinematicsUtils::momentumInLab(s, m1, m2);
  }

  G4double KinematicsUtils::sumTotalEnergies(const ParticleList &pl) {
    G4double E = 0.0;
    for(ParticleIter i = pl.begin(); i != pl.end(); ++i) {
      E += (*i)->getEnergy();
    }
    return E;
  }

  ThreeVector KinematicsUtils::sumMomenta(const ParticleList &pl) {
    ThreeVector p(0.0, 0.0, 0.0);
    for(ParticleIter i = pl.begin(); i != pl.end(); ++i) {
      p += (*i)->getMomentum();
    }
    return p;
  }

  G4double KinematicsUtils::energy(const ThreeVector &p, const G4double m) {
    return std::sqrt(p.mag2() + m*m);
  }

  G4double KinematicsUtils::invariantMass(const G4double E, const ThreeVector & p) {
    return std::sqrt(E*E - p.mag2());
  }

  G4double KinematicsUtils::gammaFromKineticEnergy(const ParticleSpecies &p, const G4double EKin) {
    G4double mass;
    if(p.theType==Composite)
      mass = ParticleTable::getTableMass(p.theA, p.theZ);
    else
      mass = ParticleTable::getTableParticleMass(p.theType);
    return (1.+EKin/mass);
  }

}
