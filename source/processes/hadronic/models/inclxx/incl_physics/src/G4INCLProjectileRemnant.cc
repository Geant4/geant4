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

/** \file G4INCLProjectileRemnant.cc
 * \brief Class for constructing a projectile-like remnant.
 *
 * \date 20 March 2012
 * \author Davide Mancusi
 */

#include "G4INCLProjectileRemnant.hh"
#include <algorithm>
#include <numeric>

namespace G4INCL {

  void ProjectileRemnant::reset() {
    deleteParticles();
    thePosition = ThreeVector();
    theMomentum = ThreeVector();
    theEnergy = 0.0;
    thePotentialEnergy = 0.0;
    theA = 0;
    theZ = 0;
    nCollisions = 0;

    for(std::map<long, Particle*>::const_iterator i=storedComponents.begin(); i!=storedComponents.end(); ++i) {
      Particle *p = new Particle(*(i->second));
      EnergyLevelMap::iterator energyIter = theInitialEnergyLevels.find(i->first);
// assert(energyIter!=theInitialEnergyLevels.end());
      const G4double energyLevel = energyIter->second;
      theInitialEnergyLevels.erase(energyIter);
      theInitialEnergyLevels[p->getID()] = energyLevel;
      addParticle(p);
    }
    if(theA>0)
      thePosition /= theA;
    setTableMass();
    INCL_DEBUG("ProjectileRemnant object was reset:" << '\n' << print());
  }

  void ProjectileRemnant::removeParticle(Particle * const p, const G4double theProjectileCorrection) {
// assert(p->isNucleon() || p->isLambda());

    INCL_DEBUG("The following Particle is about to be removed from the ProjectileRemnant:"
        << '\n' << p->print()
        << "theProjectileCorrection=" << theProjectileCorrection << '\n');
    // Update A, Z, S, momentum, and energy of the projectile remnant
    theA -= p->getA();
    theZ -= p->getZ();
    theS -= p->getS();

    ThreeVector const &oldMomentum = p->getMomentum();
    const G4double oldEnergy = p->getEnergy();
    Cluster::removeParticle(p);

#if !defined(NDEBUG) && !defined(INCLXX_IN_GEANT4_MODE)
    ThreeVector theTotalMomentum;
    G4double theTotalEnergy = 0.;
    const G4double theThreshold = 0.1;
#endif

    if(getA()>0) { // if there are any particles left
// assert((unsigned int)getA()==particles.size());

      const G4double theProjectileCorrectionPerNucleon = theProjectileCorrection / particles.size();

      // Update the kinematics of the components
      for(ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
        (*i)->setEnergy((*i)->getEnergy() + theProjectileCorrectionPerNucleon);
        (*i)->setMass((*i)->getInvariantMass());
#if !defined(NDEBUG) && !defined(INCLXX_IN_GEANT4_MODE)
        theTotalMomentum += (*i)->getMomentum();
        theTotalEnergy += (*i)->getEnergy();
#endif
      }
    }

    theMomentum -= oldMomentum;
    theEnergy -= oldEnergy - theProjectileCorrection;

// assert(std::abs((theTotalMomentum-theMomentum).mag())<theThreshold);
// assert(std::abs(theTotalEnergy-theEnergy)<theThreshold);
    INCL_DEBUG("After Particle removal, the ProjectileRemnant looks like this:"
        << '\n' << print());
  }

  ParticleList ProjectileRemnant::addDynamicalSpectators(ParticleList pL) {
    // Try as hard as possible to add back all the dynamical spectators.
    // Don't add spectators that lead to negative excitation energies, but
    // iterate over the spectators as many times as possible, until
    // absolutely sure that all of them were rejected.
    unsigned int accepted;
    unsigned long loopCounter = 0;
    const unsigned long maxLoopCounter = 10000000;
    do {
      accepted = 0;
      ParticleList toBeAdded = pL;
      for(ParticleIter p=toBeAdded.begin(), e=toBeAdded.end(); p!=e; ++p) {
        G4bool isAccepted = addDynamicalSpectator(*p);
        if(isAccepted) {
          pL.remove(*p);
          accepted++;
        }
      }
      ++loopCounter;
    } while(loopCounter<maxLoopCounter && accepted > 0); /* Loop checking, 10.07.2015, D.Mancusi */
    return pL;
  }

  ParticleList ProjectileRemnant::addAllDynamicalSpectators(ParticleList const &pL) {
    // Put all the spectators in the projectile
    ThreeVector theNewMomentum = theMomentum;
    G4double theNewEnergy = theEnergy;
    G4int theNewA = theA;
    G4int theNewZ = theZ;
    G4int theNewS = theS;
    for(ParticleIter p=pL.begin(), e=pL.end(); p!=e; ++p) {
// assert((*p)->isNucleonorLambda());
      // Add the initial (off-shell) momentum and energy to the projectile remnant
      theNewMomentum += getStoredMomentum(*p);
      theNewEnergy += (*p)->getEnergy();
      theNewA += (*p)->getA();
      theNewZ += (*p)->getZ();
      theNewS += (*p)->getS();
    }

    // Check that the excitation energy of the new projectile remnant is non-negative
    const G4double theNewMass = ParticleTable::getTableMass(theNewA,theNewZ,theNewS);
    const G4double theNewExcitationEnergy = computeExcitationEnergyWith(pL);
    const G4double theNewEffectiveMass = theNewMass + theNewExcitationEnergy;

    // If this condition is satisfied, there is no solution. Fall back on the
    // "most" method
    if(theNewEnergy<theNewEffectiveMass) {
      INCL_WARN("Could not add all the dynamical spectators back into the projectile remnant."
           << " Falling back to the \"most\" method." << '\n');
      return addMostDynamicalSpectators(pL);
    }

    // Add all the participants to the projectile remnant
    for(ParticleIter p=pL.begin(), e=pL.end(); p!=e; ++p) {
      particles.push_back(*p);
    }

    // Rescale the momentum of the projectile remnant so that sqrt(s) has the
    // correct value
    const G4double scalingFactorSquared = (theNewEnergy*theNewEnergy-theNewEffectiveMass*theNewEffectiveMass)/theNewMomentum.mag2();
    const G4double scalingFactor = std::sqrt(scalingFactorSquared);
    INCL_DEBUG("Scaling factor for the projectile-remnant momentum = " << scalingFactor << '\n');

    theA = theNewA;
    theZ = theNewZ;
    theS = theNewS;
    theMomentum = theNewMomentum * scalingFactor;
    theEnergy = theNewEnergy;

    return ParticleList();
  }

  ParticleList ProjectileRemnant::addMostDynamicalSpectators(ParticleList pL) {
    // Try as hard as possible to add back all the dynamical spectators.
    // Don't add spectators that lead to negative excitation energies. Start by
    // adding all of them, and repeatedly remove the most troublesome one until
    // the excitation energy becomes non-negative.

    // Put all the spectators in the projectile
    ThreeVector theNewMomentum = theMomentum;
    G4double theNewEnergy = theEnergy;
    G4int theNewA = theA;
    G4int theNewZ = theZ;
    G4int theNewS = theS;
    for(ParticleIter p=pL.begin(), e=pL.end(); p!=e; ++p) {
// assert((*p)->isNucleonorLambda());
      // Add the initial (off-shell) momentum and energy to the projectile remnant
      theNewMomentum += getStoredMomentum(*p);
      theNewEnergy += (*p)->getEnergy();
      theNewA += (*p)->getA();
      theNewZ += (*p)->getZ();
      theNewS += (*p)->getS();
    }

    // Check that the excitation energy of the new projectile remnant is non-negative
    const G4double theNewMass = ParticleTable::getTableMass(theNewA,theNewZ,theNewS);
    const G4double theNewInvariantMassSquared = theNewEnergy*theNewEnergy-theNewMomentum.mag2();

    G4bool positiveExcitationEnergy = false;
    if(theNewInvariantMassSquared>=0.) {
      const G4double theNewInvariantMass = std::sqrt(theNewInvariantMassSquared);
      positiveExcitationEnergy = (theNewInvariantMass-theNewMass>-1.e-5);
    }

    // Keep removing nucleons from the projectile remnant until we achieve a
    // non-negative excitation energy.
    ParticleList rejected;
    while(!positiveExcitationEnergy && !pL.empty()) { /* Loop checking, 10.07.2015, D.Mancusi */
      G4double maxExcitationEnergy = -1.E30;
      ParticleMutableIter best = pL.end();
      ThreeVector bestMomentum;
      G4double bestEnergy = -1.;
      G4int bestA = -1, bestZ = -1, bestS = 0;
      for(ParticleList::iterator p=pL.begin(), e=pL.end(); p!=e; ++p) {
        // Subtract the initial (off-shell) momentum and energy from the new
        // projectile remnant
        const ThreeVector theNewerMomentum = theNewMomentum - getStoredMomentum(*p);
        const G4double theNewerEnergy = theNewEnergy - (*p)->getEnergy();
        const G4int theNewerA = theNewA - (*p)->getA();
        const G4int theNewerZ = theNewZ - (*p)->getZ();
        const G4int theNewerS = theNewS - (*p)->getS();

        const G4double theNewerMass = ParticleTable::getTableMass(theNewerA,theNewerZ,theNewerS);
        const G4double theNewerInvariantMassSquared = theNewerEnergy*theNewerEnergy-theNewerMomentum.mag2();

        if(theNewerInvariantMassSquared>=-1.e-5) {
          const G4double theNewerInvariantMass = std::sqrt(std::max(0.,theNewerInvariantMassSquared));
          const G4double theNewerExcitationEnergy = ((theNewerA>1) ? theNewerInvariantMass-theNewerMass : 0.);
          // Pick the nucleon that maximises the excitation energy of the
          // ProjectileRemnant
          if(theNewerExcitationEnergy>maxExcitationEnergy) {
            best = p;
            maxExcitationEnergy = theNewerExcitationEnergy;
            bestMomentum = theNewerMomentum;
            bestEnergy = theNewerEnergy;
            bestA = theNewerA;
            bestZ = theNewerZ;
            bestS = theNewerS;
          }
        }
      }

      // If we couldn't even calculate the excitation energy, fail miserably
      if(best==pL.end())
        return pL;

      rejected.push_back(*best);
      pL.erase(best);
      theNewMomentum = bestMomentum;
      theNewEnergy = bestEnergy;
      theNewA = bestA;
      theNewZ = bestZ;
      theNewS = bestS;

      if(maxExcitationEnergy>0.) {
        // Stop here
        positiveExcitationEnergy = true;
      }
    }

    // Add the accepted participants to the projectile remnant
    for(ParticleIter p=pL.begin(), e=pL.end(); p!=e; ++p) {
      particles.push_back(*p);
    }
    theA = theNewA;
    theZ = theNewZ;
    theS = theNewS;
    theMomentum = theNewMomentum;
    theEnergy = theNewEnergy;

    return rejected;
  }

  G4bool ProjectileRemnant::addDynamicalSpectator(Particle * const p) {
// assert(p->isNucleon());

    // Add the initial (off-shell) momentum and energy to the projectile remnant
    ThreeVector const &oldMomentum = getStoredMomentum(p);
    const ThreeVector theNewMomentum = theMomentum + oldMomentum;
    const G4double oldEnergy = p->getEnergy();
    const G4double theNewEnergy = theEnergy + oldEnergy;

    // Check that the excitation energy of the new projectile remnant is non-negative
    const G4double theNewMass = ParticleTable::getTableMass(theA+p->getA(),theZ+p->getZ(),theS+p->getS());
    const G4double theNewInvariantMassSquared = theNewEnergy*theNewEnergy-theNewMomentum.mag2();

    if(theNewInvariantMassSquared<0.)
      return false;

    const G4double theNewInvariantMass = std::sqrt(theNewInvariantMassSquared);

    if(theNewInvariantMass-theNewMass<-1.e-5)
      return false; // negative excitation energy here

    // Add the spectator to the projectile remnant
    theA += p->getA();
    theZ += p->getZ();
    theMomentum = theNewMomentum;
    theEnergy = theNewEnergy;
    particles.push_back(p);
    return true;
  }

  G4double ProjectileRemnant::computeExcitationEnergyExcept(const long exceptID) const {
    const EnergyLevels theEnergyLevels = getPresentEnergyLevelsExcept(exceptID);
    return computeExcitationEnergy(theEnergyLevels);
  }

  G4double ProjectileRemnant::computeExcitationEnergyWith(const ParticleList &pL) const {
    const EnergyLevels theEnergyLevels = getPresentEnergyLevelsWith(pL);
    return computeExcitationEnergy(theEnergyLevels);
  }

  G4double ProjectileRemnant::computeExcitationEnergy(const EnergyLevels &levels) const {
    // The ground-state energy is the sum of the A smallest initial projectile
    // energies.
    // For the last nucleon, return 0 so that the algorithm will just put it on
    // shell.
    const std::size_t theNewA = levels.size();
// assert(theNewA>0);
    if(theNewA==1)
      return 0.;

    const G4double groundState = theGroundStateEnergies.at(theNewA-1);

    // Compute the sum of the presently occupied energy levels
    const G4double excitedState = std::accumulate(
        levels.cbegin(),
        levels.cend(),
        0.);

    return excitedState-groundState;
  }

  ProjectileRemnant::EnergyLevels ProjectileRemnant::getPresentEnergyLevelsExcept(const long exceptID) const {
    EnergyLevels theEnergyLevels;
    for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
      if((*p)->getID()!=exceptID) {
        EnergyLevelMap::const_iterator i = theInitialEnergyLevels.find((*p)->getID());
// assert(i!=theInitialEnergyLevels.end());
        theEnergyLevels.push_back(i->second);
      }
    }
// assert(theEnergyLevels.size()==particles.size()-1);
    return theEnergyLevels;
  }

  ProjectileRemnant::EnergyLevels ProjectileRemnant::getPresentEnergyLevelsWith(const ParticleList &pL) const {
    EnergyLevels theEnergyLevels;
    for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
      EnergyLevelMap::const_iterator i = theInitialEnergyLevels.find((*p)->getID());
// assert(i!=theInitialEnergyLevels.end());
      theEnergyLevels.push_back(i->second);
    }
    for(ParticleIter p=pL.begin(), e=pL.end(); p!=e; ++p) {
      EnergyLevelMap::const_iterator i = theInitialEnergyLevels.find((*p)->getID());
// assert(i!=theInitialEnergyLevels.end());
      theEnergyLevels.push_back(i->second);
    }

// assert(theEnergyLevels.size()==particles.size()+pL.size());
    return theEnergyLevels;
  }

}

