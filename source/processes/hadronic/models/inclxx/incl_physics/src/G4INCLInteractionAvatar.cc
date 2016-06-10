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

/* \file G4INCLInteractionAvatar.cc
 * \brief Virtual class for interaction avatars.
 *
 * This class is inherited by decay and collision avatars. The goal is to
 * provide a uniform treatment of common physics, such as Pauli blocking,
 * enforcement of energy conservation, etc.
 *
 *  \date Mar 1st, 2011
 * \author Davide Mancusi
 */

#include "G4INCLInteractionAvatar.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLPauliBlocking.hh"
#include "G4INCLRootFinder.hh"
#include "G4INCLLogger.hh"
#include "G4INCLConfigEnums.hh"
// #include <cassert>

namespace G4INCL {

  const G4double InteractionAvatar::locEAccuracy = 1.E-4;
  const G4int InteractionAvatar::maxIterLocE = 50;
  G4ThreadLocal Particle *InteractionAvatar::backupParticle1 = NULL;
  G4ThreadLocal Particle *InteractionAvatar::backupParticle2 = NULL;

  InteractionAvatar::InteractionAvatar(G4double time, G4INCL::Nucleus *n, G4INCL::Particle *p1)
    : IAvatar(time), theNucleus(n),
    particle1(p1), particle2(NULL),
    isPiN(false),
    violationEFunctor(NULL)
  {
  }

  InteractionAvatar::InteractionAvatar(G4double time, G4INCL::Nucleus *n, G4INCL::Particle *p1,
      G4INCL::Particle *p2)
    : IAvatar(time), theNucleus(n),
    particle1(p1), particle2(p2),
    isPiN((p1->isPion() && p2->isNucleon()) || (p2->isPion() && p1->isNucleon())),
    violationEFunctor(NULL)
  {
  }

  InteractionAvatar::~InteractionAvatar() {
  }

  void InteractionAvatar::deleteBackupParticles() {
    delete backupParticle1;
    if(backupParticle2)
      delete backupParticle2;
    backupParticle1 = NULL;
    backupParticle2 = NULL;
  }

  void InteractionAvatar::preInteractionBlocking() {
    if(backupParticle1)
      (*backupParticle1) = (*particle1);
    else
      backupParticle1 = new Particle(*particle1);

    if(particle2) {
      if(backupParticle2)
        (*backupParticle2) = (*particle2);
      else
        backupParticle2 = new Particle(*particle2);

      oldTotalEnergy = particle1->getEnergy() + particle2->getEnergy()
        - particle1->getPotentialEnergy() - particle2->getPotentialEnergy();
      oldXSec = CrossSections::total(particle1, particle2);
    } else {
      oldTotalEnergy = particle1->getEnergy() - particle1->getPotentialEnergy();
    }
  }

  void InteractionAvatar::preInteractionLocalEnergy(Particle * const p) {
    if(!theNucleus || p->isPion()) return; // Local energy does not make any sense without a nucleus

    if(shouldUseLocalEnergy())
      KinematicsUtils::transformToLocalEnergyFrame(theNucleus, p);
  }

  void InteractionAvatar::preInteraction() {
    preInteractionBlocking();

    preInteractionLocalEnergy(particle1);

    if(particle2) {
      preInteractionLocalEnergy(particle2);
      if(!isPiN) {
        boostVector = KinematicsUtils::makeBoostVector(particle1, particle2);
        particle2->boost(boostVector);
      }
    } else {
      boostVector = particle1->getMomentum()/particle1->getEnergy();
    }
    if(!isPiN)
      particle1->boost(boostVector);
  }

  G4bool InteractionAvatar::bringParticleInside(Particle * const p) {
    if(!theNucleus)
      return false;

    ThreeVector pos = p->getPosition();
    p->rpCorrelate();
    G4double pos2 = pos.mag2();
    const G4double r = theNucleus->getSurfaceRadius(p);
    short iterations=0;
    const short maxIterations=50;

    if(pos2 < r*r) return true;

    while( pos2 >= r*r && iterations<maxIterations )
    {
      pos *= std::sqrt(r*r*0.99/pos2);
      pos2 = pos.mag2();
      iterations++;
    }
    if( iterations < maxIterations)
    {
      INCL_DEBUG("Particle position vector length was : " << p->getPosition().mag() << ", rescaled to: " << pos.mag() << std::endl);
      p->setPosition(pos);
      return true;
    }
    else
      return false;
  }

  FinalState *InteractionAvatar::postInteraction(FinalState *fs) {
    INCL_DEBUG("postInteraction: final state: " << std::endl << fs->print() << std::endl);
    ParticleList modified = fs->getModifiedParticles();
    ParticleList modifiedAndCreated = modified;
    ParticleList created = fs->getCreatedParticles();
    modifiedAndCreated.insert(modifiedAndCreated.end(), created.begin(), created.end());

    if(!isPiN) {
      // Boost back to lab
      for(ParticleIter i=modifiedAndCreated.begin(), e=modifiedAndCreated.end(); i!=e; ++i )
        (*i)->boost(-boostVector);
    }

    // If there is no Nucleus, just return
    if(!theNucleus) return fs;

    // Mark pions that have been created outside their well (we will force them
    // to be emitted later).
    for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i )
      if((*i)->isPion() && (*i)->getPosition().mag() > theNucleus->getSurfaceRadius(*i)) {
        (*i)->makeParticipant();
        (*i)->setOutOfWell();
        fs->addOutgoingParticle(*i);
        INCL_DEBUG("Pion was created outside its potential well." << std::endl
            << (*i)->print());
      }

    // Try to enforce energy conservation
    fs->setTotalEnergyBeforeInteraction(oldTotalEnergy);
    G4bool success = true;
    if(!isPiN || shouldUseLocalEnergy())
      success = enforceEnergyConservation(fs);
    if(!success) {
      INCL_DEBUG("Enforcing energy conservation: failed!" << std::endl);

      // Restore the state of the initial particles
      restoreParticles();

      // Delete newly created particles
      for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i )
        delete *i;

      FinalState *fsBlocked = new FinalState;
      delete fs;
      fsBlocked->makeNoEnergyConservation();
      fsBlocked->setTotalEnergyBeforeInteraction(0.0);

      return fsBlocked; // Interaction is blocked. Return an empty final state.
    }
    INCL_DEBUG("Enforcing energy conservation: success!" << std::endl);

    INCL_DEBUG("postInteraction after energy conservation: final state: " << std::endl << fs->print() << std::endl);

    // Check that outgoing delta resonances can decay to pi-N
    for(ParticleIter i=modified.begin(), e=modified.end(); i!=e; ++i )
      if((*i)->isDelta() &&
          (*i)->getMass() < ParticleTable::effectiveDeltaDecayThreshold) {
        INCL_DEBUG("Mass of the produced delta below decay threshold; forbidding collision. deltaMass=" <<
            (*i)->getMass() << std::endl);

        // Restore the state of the initial particles
        restoreParticles();

        // Delete newly created particles
        for(ParticleIter j=created.begin(), end=created.end(); j!=end; ++j )
          delete *j;

        FinalState *fsBlocked = new FinalState;
        delete fs;
        fsBlocked->makeNoEnergyConservation();
        fsBlocked->setTotalEnergyBeforeInteraction(0.0);

        return fsBlocked; // Interaction is blocked. Return an empty final state.
      }

    INCL_DEBUG("Random seeds before Pauli blocking: " << Random::getSeeds() << std::endl);
    // Test Pauli blocking
    G4bool isBlocked = Pauli::isBlocked(modifiedAndCreated, theNucleus);

    if(isBlocked) {
      INCL_DEBUG("Pauli: Blocked!" << std::endl);

      // Restore the state of the initial particles
      restoreParticles();

      // Delete newly created particles
      for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i )
        delete *i;

      FinalState *fsBlocked = new FinalState;
      delete fs;
      fsBlocked->makePauliBlocked();
      fsBlocked->setTotalEnergyBeforeInteraction(0.0);

      return fsBlocked; // Interaction is blocked. Return an empty final state.
    }
    INCL_DEBUG("Pauli: Allowed!" << std::endl);

    // Test CDPP blocking
    G4bool isCDPPBlocked = Pauli::isCDPPBlocked(created, theNucleus);

    if(isCDPPBlocked) {
      INCL_DEBUG("CDPP: Blocked!" << std::endl);

      // Restore the state of the initial particles
      restoreParticles();

      // Delete newly created particles
      for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i )
        delete *i;

      FinalState *fsBlocked = new FinalState;
      delete fs;
      fsBlocked->makePauliBlocked();
      fsBlocked->setTotalEnergyBeforeInteraction(0.0);

      return fsBlocked; // Interaction is blocked. Return an empty final state.
    }
    INCL_DEBUG("CDPP: Allowed!" << std::endl);

    // If all went well, try to bring particles inside the nucleus...
    for(ParticleIter i=modifiedAndCreated.begin(), e=modifiedAndCreated.end(); i!=e; ++i )
    {
      // ...except for pions beyond their surface radius.
      if((*i)->isOutOfWell()) continue;

      const G4bool successBringParticlesInside = bringParticleInside(*i);
      if( !successBringParticlesInside ) {
        INCL_ERROR("Failed to bring particle inside the nucleus!" << std::endl);
      }
    }

    // Collision accepted!
    for(ParticleIter i=modifiedAndCreated.begin(), e=modifiedAndCreated.end(); i!=e; ++i ) {
      if(!(*i)->isOutOfWell()) {
        // Decide if the particle should be made into a spectator
        // (Back to spectator)
        G4bool goesBackToSpectator = false;
        if((*i)->isNucleon() && theNucleus->getStore()->getConfig()->getBackToSpectator()) {
          G4double threshold = (*i)->getPotentialEnergy();
          if((*i)->getType()==Proton)
            threshold += Math::twoThirds*theNucleus->getTransmissionBarrier(*i);
          if((*i)->getKineticEnergy() < threshold)
            goesBackToSpectator = true;
        }

        // Thaw the particle propagation
        (*i)->thawPropagation();

        // Increment or decrement the participant counters
        if(goesBackToSpectator) {
          INCL_DEBUG("The following particle goes back to spectator:" << std::endl
              << (*i)->print() << std::endl);
          if(!(*i)->isTargetSpectator()) {
            theNucleus->getStore()->getBook().decrementCascading();
          }
          (*i)->makeTargetSpectator();
        } else {
          if((*i)->isTargetSpectator()) {
            theNucleus->getStore()->getBook().incrementCascading();
          }
          (*i)->makeParticipant();
        }
      }
    }
    ParticleList destroyed = fs->getDestroyedParticles();
    for(ParticleIter i=destroyed.begin(), e=destroyed.end(); i!=e; ++i )
      if(!(*i)->isTargetSpectator())
        theNucleus->getStore()->getBook().decrementCascading();

    return fs;
  }

  void InteractionAvatar::restoreParticles() const {
    (*particle1) = (*backupParticle1);
    if(particle2)
      (*particle2) = (*backupParticle2);
  }

  G4bool InteractionAvatar::enforceEnergyConservation(FinalState * const fs) {
    // Set up the violationE calculation
    ParticleList modified = fs->getModifiedParticles();
    const G4bool manyBodyFinalState = (modified.size() + fs->getCreatedParticles().size() > 1);
    if(manyBodyFinalState)
      violationEFunctor = new ViolationEMomentumFunctor(theNucleus, fs, &boostVector, shouldUseLocalEnergy());
    else {
      Particle const * const p = modified.front();
      // The following condition is necessary for the functor to work
      // correctly. A similar condition exists in INCL4.6.
      if(p->getMass() < ParticleTable::effectiveDeltaDecayThreshold)
        return false;
      violationEFunctor = new ViolationEEnergyFunctor(theNucleus, fs, shouldUseLocalEnergy());
    }

    // Apply the root-finding algorithm
    const RootFinder::Solution theSolution = RootFinder::solve(violationEFunctor, 1.0);
    if(theSolution.success) { // Apply the solution
      (*violationEFunctor)(theSolution.x);
    } else if(theNucleus){
      INCL_DEBUG("Couldn't enforce energy conservation after an interaction, root-finding algorithm failed." << std::endl);
      theNucleus->getStore()->getBook().incrementEnergyViolationInteraction();
    }
    delete violationEFunctor;
    violationEFunctor = NULL;
    return theSolution.success;
  }

  /* ***                                                      ***
   * *** InteractionAvatar::ViolationEMomentumFunctor methods ***
   * ***                                                      ***/

  InteractionAvatar::ViolationEMomentumFunctor::ViolationEMomentumFunctor(Nucleus * const nucleus, FinalState const * const finalState, ThreeVector const * const boost, const G4bool localE) :
    RootFunctor(0., 1E6),
    initialEnergy(finalState->getTotalEnergyBeforeInteraction()),
    theNucleus(nucleus),
    boostVector(boost),
    shouldUseLocalEnergy(localE)
  {
    // Set up the finalParticles list
    finalParticles = finalState->getModifiedParticles();
    ParticleList const &created = finalState->getCreatedParticles();
    finalParticles.insert(finalParticles.end(), created.begin(), created.end());

    // Store the particle momenta (necessary for the calls to
    // scaleParticleMomenta() to work)
    particleMomenta.clear();
    for(ParticleIter i=finalParticles.begin(), e=finalParticles.end(); i!=e; ++i) {
      (*i)->boost(*boostVector);
      particleMomenta.push_back((*i)->getMomentum());
    }
  }

  G4double InteractionAvatar::ViolationEMomentumFunctor::operator()(const G4double alpha) const {
    scaleParticleMomenta(alpha);

    G4double deltaE = 0.0;
    for(ParticleIter i=finalParticles.begin(), e=finalParticles.end(); i!=e; ++i)
      deltaE += (*i)->getEnergy() - (*i)->getPotentialEnergy();
    deltaE -= initialEnergy;
    return deltaE;
  }

  void InteractionAvatar::ViolationEMomentumFunctor::scaleParticleMomenta(const G4double alpha) const {

    std::list<ThreeVector>::const_iterator iP = particleMomenta.begin();
    for(ParticleIter i=finalParticles.begin(), e=finalParticles.end(); i!=e; ++i, ++iP) {
      (*i)->setMomentum((*iP)*alpha);
      (*i)->adjustEnergyFromMomentum();
      (*i)->rpCorrelate();
      (*i)->boost(-(*boostVector));
      if(theNucleus)
        theNucleus->updatePotentialEnergy(*i);
      else
        (*i)->setPotentialEnergy(0.);

      if(shouldUseLocalEnergy && !(*i)->isPion()) { // This translates AECSVT's loops 1, 3 and 4
// assert(theNucleus); // Local energy without a nucleus doesn't make sense
        const G4double energy = (*i)->getEnergy(); // Store the energy of the particle
        G4double locE = KinematicsUtils::getLocalEnergy(theNucleus, *i); // Initial value of local energy
        G4double locEOld;
        G4double deltaLocE = InteractionAvatar::locEAccuracy + 1E3;
        for(G4int iterLocE=0;
            deltaLocE>InteractionAvatar::locEAccuracy && iterLocE<InteractionAvatar::maxIterLocE;
            ++iterLocE) {
          locEOld = locE;
          (*i)->setEnergy(energy + locE); // Update the energy of the particle...
          (*i)->adjustMomentumFromEnergy();
          theNucleus->updatePotentialEnergy(*i); // ...update its potential energy...
          locE = KinematicsUtils::getLocalEnergy(theNucleus, *i); // ...and recompute locE.
          deltaLocE = std::abs(locE-locEOld);
        }
      }
    }
  }

  void InteractionAvatar::ViolationEMomentumFunctor::cleanUp(const G4bool success) const {
    if(!success)
      scaleParticleMomenta(1.);
  }

  /* ***                                                    ***
   * *** InteractionAvatar::ViolationEEnergyFunctor methods ***
   * ***                                                    ***/

  InteractionAvatar::ViolationEEnergyFunctor::ViolationEEnergyFunctor(Nucleus * const nucleus, FinalState const * const finalState, const G4bool localE) :
    RootFunctor(0., 1E6),
    initialEnergy(finalState->getTotalEnergyBeforeInteraction()),
    theNucleus(nucleus),
    theParticle(finalState->getModifiedParticles().front()),
    theEnergy(theParticle->getEnergy()),
    theMomentum(theParticle->getMomentum()),
    energyThreshold(KinematicsUtils::energy(theMomentum,ParticleTable::effectiveDeltaDecayThreshold)),
    shouldUseLocalEnergy(localE)
  {
// assert(finalState->getModifiedParticles().size()==1);
// assert(theParticle->isDelta());
  }

  G4double InteractionAvatar::ViolationEEnergyFunctor::operator()(const G4double alpha) const {
    setParticleEnergy(alpha);
    return theParticle->getEnergy() - theParticle->getPotentialEnergy() - initialEnergy;
  }

  void InteractionAvatar::ViolationEEnergyFunctor::setParticleEnergy(const G4double alpha) const {

    G4double locE;
    if(shouldUseLocalEnergy) {
// assert(theNucleus); // Local energy without a nucleus doesn't make sense
      locE = KinematicsUtils::getLocalEnergy(theNucleus, theParticle); // Initial value of local energy
    } else
      locE = 0.;
    G4double locEOld;
    G4double deltaLocE = InteractionAvatar::locEAccuracy + 1E3;
    for(G4int iterLocE=0;
        deltaLocE>InteractionAvatar::locEAccuracy && iterLocE<InteractionAvatar::maxIterLocE;
        ++iterLocE) {
      locEOld = locE;
      const G4double particleEnergy = energyThreshold + alpha*(theEnergy-energyThreshold);
      const G4double theMass = std::sqrt(std::pow(particleEnergy,2.)-theMomentum.mag2());
      theParticle->setMass(theMass);
      theParticle->setEnergy(particleEnergy + locE); // Update the energy of the particle...
      theParticle->adjustMomentumFromEnergy();
      if(theNucleus) {
        theNucleus->updatePotentialEnergy(theParticle); // ...update its potential energy...
        if(shouldUseLocalEnergy)
          locE = KinematicsUtils::getLocalEnergy(theNucleus, theParticle); // ...and recompute locE.
        else
          locE = 0.;
      } else
        locE = 0.;
      deltaLocE = std::abs(locE-locEOld);
    }

  }

  void InteractionAvatar::ViolationEEnergyFunctor::cleanUp(const G4bool success) const {
    if(!success)
      setParticleEnergy(1.);
  }

}
