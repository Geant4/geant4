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

#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"

namespace G4INCL {

  namespace KinematicsUtils {

  void transformToLocalEnergyFrame(Nucleus const * const n, Particle * const p) {
// assert(!p->isMeson()); // No local energy for mesons
    const G4double localEnergy = getLocalEnergy(n, p);
    const G4double localTotalEnergy = p->getEnergy() - localEnergy;
    p->setEnergy(localTotalEnergy);
    p->adjustMomentumFromEnergy();
  }

  G4double getLocalEnergy(Nucleus const * const n, Particle * const p) {
// assert(!p->isMeson()); // No local energy for mesons
    
    G4double vloc = 0.0;
    const G4double r = p->getPosition().mag();
    const G4double mass = p->getMass();

    // Local energy is constant outside the surface
    if(r > n->getUniverseRadius()) {
      INCL_WARN("Tried to evaluate local energy for a particle outside the maximum radius."
            << '\n' << p->print() << '\n'
            << "Maximum radius = " << n->getDensity()->getMaximumRadius() << '\n'
            << "Universe radius = " << n->getUniverseRadius() << '\n');
      return 0.0;
    }

    G4double pfl0 = 0.0;
    const ParticleType t = p->getType();
    const G4double kinE = p->getKineticEnergy();
    if(kinE <= n->getPotential()->getFermiEnergy(t)) {
      pfl0 = n->getPotential()->getFermiMomentum(p);
    } else {
      const G4double tf0 = p->getPotentialEnergy() - n->getPotential()->getSeparationEnergy(p);
      if(tf0<0.0) return 0.0;
      pfl0 = std::sqrt(tf0*(tf0 + 2.0*mass));
    }
    const G4double pReflection = p->getReflectionMomentum()/pfl0;
    const G4double reflectionRadius = n->getDensity()->getMaxRFromP(p->getType(), pReflection);
    const G4double pNominal = p->getMomentum().mag()/pfl0;
    const G4double nominalReflectionRadius = n->getDensity()->getMaxRFromP(p->getType(), pNominal);
    const G4double pl = pfl0*n->getDensity()->getMinPFromR(t, r*nominalReflectionRadius/reflectionRadius);
    vloc = std::sqrt(pl*pl + mass*mass) - mass;

    return vloc;
  }

  ThreeVector makeBoostVector(Particle const * const p1, Particle const * const p2){
    const G4double totalEnergy = p1->getEnergy() + p2->getEnergy();
    return ((p1->getMomentum() + p2->getMomentum())/totalEnergy);
  }

  G4double totalEnergyInCM(Particle const * const p1, Particle const * const p2){
    return std::sqrt(squareTotalEnergyInCM(p1,p2));
  }

  G4double squareTotalEnergyInCM(Particle const * const p1, Particle const * const p2) {
    G4double beta2 = makeBoostVector(p1, p2).mag2();
    if(beta2 > 1.0) {
      INCL_ERROR("squareTotalEnergyInCM: beta2 == " << beta2 << " > 1.0" << '\n');
      beta2 = 0.0;
    }
    return (1.0 - beta2)*std::pow(p1->getEnergy() + p2->getEnergy(), 2);
  }

  G4double momentumInCM(Particle const * const p1, Particle const * const p2) {
    const G4double m1sq = std::pow(p1->getMass(),2);
    const G4double m2sq = std::pow(p2->getMass(),2);
    const G4double z = p1->getEnergy()*p2->getEnergy() - p1->getMomentum().dot(p2->getMomentum());
    G4double pcm2 = (z*z-m1sq*m2sq)/(2*z+m1sq+m2sq);
    if(pcm2 < 0.0) {
      INCL_ERROR("momentumInCM: pcm2 == " << pcm2 << " < 0.0" << '\n');
      pcm2 = 0.0;
    }
    return std::sqrt(pcm2);
  }

  G4double momentumInCM(const G4double E, const G4double M1, const G4double M2) {
    return 0.5*std::sqrt((E*E - std::pow(M1 + M2, 2))
			 *(E*E - std::pow(M1 - M2, 2)))/E;
  }

  G4double momentumInLab(const G4double s, const G4double m1, const G4double m2) {
    const G4double m1sq = m1*m1;
    const G4double m2sq = m2*m2;
    G4double plab2 = (s*s-2*s*(m1sq+m2sq)+(m1sq-m2sq)*(m1sq-m2sq))/(4*m2sq);
    if(plab2 < 0.0) {
      INCL_ERROR("momentumInLab: plab2 == " << plab2 << " < 0.0; m1sq == " << m1sq << "; m2sq == " << m2sq << "; s == " << s << '\n');
      plab2 = 0.0;
    }
    return std::sqrt(plab2);
  }

  G4double momentumInLab(Particle const * const p1, Particle const * const p2) {
    const G4double m1 = p1->getMass();
    const G4double m2 = p2->getMass();
    const G4double s = squareTotalEnergyInCM(p1, p2);
    return momentumInLab(s, m1, m2);
  }

  G4double sumTotalEnergies(const ParticleList &pl) {
    G4double E = 0.0;
    for(ParticleIter i=pl.begin(), e=pl.end(); i!=e; ++i) {
      E += (*i)->getEnergy();
    }
    return E;
  }

  ThreeVector sumMomenta(const ParticleList &pl) {
    ThreeVector p(0.0, 0.0, 0.0);
    for(ParticleIter i=pl.begin(), e=pl.end(); i!=e; ++i) {
      p += (*i)->getMomentum();
    }
    return p;
  }

  G4double energy(const ThreeVector &p, const G4double m) {
    return std::sqrt(p.mag2() + m*m);
  }

  G4double invariantMass(const G4double E, const ThreeVector & p) {
    return std::sqrt(squareInvariantMass(E, p));
  }

  G4double squareInvariantMass(const G4double E, const ThreeVector & p) {
    return E*E - p.mag2();
  }

  G4double gammaFromKineticEnergy(const ParticleSpecies &p, const G4double EKin) {
    G4double mass;
    if(p.theType==Composite)
      mass = ParticleTable::getTableMass(p.theA, p.theZ);
    else
      mass = ParticleTable::getTableParticleMass(p.theType);
    return (1.+EKin/mass);
  }

  }

}
