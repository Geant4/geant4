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

  G4int shuffleComponentsHelper(G4int range) {
    return (G4int)(Random::shoot1()*range);
  }

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
    thePosition /= theA;
    setTableMass();
    DEBUG("ProjectileRemnant object was reset:" << std::endl << print());
  }

  void ProjectileRemnant::removeParticle(Particle * const p, const G4double theProjectileCorrection) {
// assert(p->isNucleon());

    DEBUG("The following Particle is about to be removed from the ProjectileRemnant:"
        << std::endl << p->print()
        << "theProjectileCorrection=" << theProjectileCorrection << std::endl);
    // Update A, Z, momentum and energy of the projectile remnant
    theA -= p->getA();
    theZ -= p->getZ();

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
      for(ParticleIter i=particles.begin(); i!=particles.end(); ++i) {
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
    DEBUG("After Particle removal, the ProjectileRemnant looks like this:"
        << std::endl << print());
  }

  ParticleList ProjectileRemnant::addDynamicalSpectators(ParticleList pL) {
    // Try as hard as possible to add back all the dynamical spectators.
    // Don't add spectators that lead to negative excitation energies, but
    // iterate over the spectators as many times as possible, until
    // absolutely sure that all of them were rejected.
    unsigned int accepted;
    do {
      accepted = 0;
      ParticleList toBeAdded = pL;
      for(ParticleIter p=toBeAdded.begin(); p!=toBeAdded.end(); ++p) {
        G4bool isAccepted = addDynamicalSpectator(*p);
        if(isAccepted) {
          pL.remove(*p);
          accepted++;
        }
      }
    } while(accepted > 0);
    return pL;
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
    for(ParticleIter p=pL.begin(); p!=pL.end(); ++p) {
// assert((*p)->isNucleon());
      // Add the initial (off-shell) momentum and energy to the projectile remnant
      theNewMomentum += getStoredMomentum(*p);
      theNewEnergy += (*p)->getEnergy();
      theNewA += (*p)->getA();
      theNewZ += (*p)->getZ();
    }

    // Check that the excitation energy of the new projectile remnant is non-negative
    const G4double theNewMass = ParticleTable::getTableMass(theNewA,theNewZ);
    const G4double theNewInvariantMassSquared = theNewEnergy*theNewEnergy-theNewMomentum.mag2();

    G4bool positiveExcitationEnergy = false;
    if(theNewInvariantMassSquared>=0.) {
      const G4double theNewInvariantMass = std::sqrt(theNewInvariantMassSquared);
      positiveExcitationEnergy = (theNewInvariantMass-theNewMass>-1.e-5);
    }

    // Keep removing nucleons from the projectile remnant until we achieve a
    // non-negative excitation energy.
    ParticleList rejected;
    while(!positiveExcitationEnergy && !pL.empty()) {
      G4double maxExcitationEnergy = -1.E30;
      ParticleList::iterator best = pL.end();
      ThreeVector bestMomentum;
      G4double bestEnergy = -1.;
      G4int bestA = -1, bestZ = -1;
      for(ParticleList::iterator p=pL.begin(); p!=pL.end(); ++p) {
        // Subtract the initial (off-shell) momentum and energy from the new
        // projectile remnant
        const ThreeVector theNewerMomentum = theNewMomentum - getStoredMomentum(*p);
        const G4double theNewerEnergy = theNewEnergy - (*p)->getEnergy();
        const G4int theNewerA = theNewA - (*p)->getA();
        const G4int theNewerZ = theNewZ - (*p)->getZ();

        const G4double theNewerMass = ParticleTable::getTableMass(theNewerA,theNewerZ);
        const G4double theNewerInvariantMassSquared = theNewerEnergy*theNewerEnergy-theNewerMomentum.mag2();

        if(theNewerInvariantMassSquared>=-1.e-5) {
          const G4double theNewerInvariantMass = std::sqrt(std::max(0.,theNewerInvariantMassSquared));
          const G4double theNewerExcitationEnergy = theNewerInvariantMass-theNewerMass;
          // Pick the nucleon that maximises the excitation energy of the
          // ProjectileRemnant
          if(theNewerExcitationEnergy>maxExcitationEnergy) {
            best = p;
            maxExcitationEnergy = theNewerExcitationEnergy;
            bestMomentum = theNewerMomentum;
            bestEnergy = theNewerEnergy;
            bestA = theNewerA;
            bestZ = theNewerZ;
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

      if(maxExcitationEnergy>0.) {
        // Stop here
        positiveExcitationEnergy = true;
      }
    }

    // Add the accepted participants to the projectile remnant
    for(ParticleIter p=pL.begin(); p!=pL.end(); ++p) {
      particles.push_back(*p);
    }
    theA = theNewA;
    theZ = theNewZ;
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
    const G4double theNewMass = ParticleTable::getTableMass(theA+p->getA(),theZ+p->getZ());
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

  G4double ProjectileRemnant::computeExcitationEnergy(const long exceptID) const {
    // The ground-state energy is the sum of the A smallest initial projectile
    // energies.
    // For the last nucleon, return 0 so that the algorithm will just put it on
    // shell.
    if(theA==1)
      return 0.;

    const G4double groundState = theGroundStateEnergies.at(theA-2);

    // Compute the sum of the presently occupied energy levels
    const EnergyLevels theEnergyLevels = getPresentEnergyLevels(exceptID);
    const G4double excitedState = std::accumulate(
        theEnergyLevels.begin(),
        theEnergyLevels.end(),
        0.);

    return excitedState-groundState;
  }

}

