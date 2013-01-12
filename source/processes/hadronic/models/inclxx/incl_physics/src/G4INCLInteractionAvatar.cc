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

  InteractionAvatar::InteractionAvatar(G4double time, G4INCL::Nucleus *n, G4INCL::Particle *p1)
    : IAvatar(time), theNucleus(n),
    particle1(p1), particle2(NULL), isPiN(false)
  {
  }

  InteractionAvatar::InteractionAvatar(G4double time, G4INCL::Nucleus *n, G4INCL::Particle *p1,
      G4INCL::Particle *p2)
    : IAvatar(time), theNucleus(n),
    particle1(p1), particle2(p2),
    isPiN((p1->isPion() && p2->isNucleon()) || (p2->isPion() && p1->isNucleon()))
  {
  }

  InteractionAvatar::~InteractionAvatar() {
  }

  void InteractionAvatar::preInteractionBlocking() {
    oldParticle1Type = particle1->getType();
    oldParticle1Energy = particle1->getEnergy();
    oldParticle1Potential = particle1->getPotentialEnergy();
    oldParticle1Momentum = particle1->getMomentum();
    oldParticle1Position = particle1->getPosition();
    oldParticle1Mass = particle1->getMass();
    oldParticle1Helicity = particle1->getHelicity();

    if(particle2) {
      oldParticle2Type = particle2->getType();
      oldParticle2Energy = particle2->getEnergy();
      oldParticle2Potential = particle2->getPotentialEnergy();
      oldParticle2Momentum = particle2->getMomentum();
      oldParticle2Position = particle2->getPosition();
      oldParticle2Mass = particle2->getMass();
      oldParticle2Helicity = particle2->getHelicity();
      oldTotalEnergy = oldParticle1Energy + oldParticle2Energy
        - particle1->getPotentialEnergy() - particle2->getPotentialEnergy();
      oldXSec = CrossSections::total(particle1, particle2);
    } else {
      oldTotalEnergy = oldParticle1Energy - particle1->getPotentialEnergy();
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
    ThreeVector pos = p->getPosition();
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
      DEBUG("Particle position vector length was : " << p->getPosition().mag() << ", rescaled to: " << pos.mag() << std::endl);
      p->setPosition(pos);
      return true;
    }
    else
      return false;
  }

  FinalState *InteractionAvatar::postInteraction(FinalState *fs) {
    ParticleList modified = fs->getModifiedParticles();
    ParticleList modifiedAndCreated = modified;
    ParticleList created = fs->getCreatedParticles();
    modifiedAndCreated.insert(modifiedAndCreated.end(), created.begin(), created.end());

    if(!isPiN) {
      // Boost back to lab
      for( ParticleIter i = modifiedAndCreated.begin(); i != modifiedAndCreated.end(); ++i )
        (*i)->boost(-boostVector);
    }

    // If there is no Nucleus, just return
    if(!theNucleus) return fs;

    // Mark pions that have been created outside their well (we will force them
    // to be emitted later).
    for( ParticleIter i = created.begin(); i != created.end(); ++i )
      if((*i)->isPion() && (*i)->getPosition().mag() > theNucleus->getSurfaceRadius(*i)) {
        (*i)->makeParticipant();
        (*i)->setOutOfWell();
        fs->addOutgoingParticle(*i);
        DEBUG("Pion was created outside its potential well." << std::endl
            << (*i)->print());
      }

    // Try to enforce energy conservation
    fs->setTotalEnergyBeforeInteraction(oldTotalEnergy);
    G4bool success = true;
    if(!isPiN || shouldUseLocalEnergy())
      success = enforceEnergyConservation(fs);
    if(!success) {
      DEBUG("Enforcing energy conservation: failed!" << std::endl);

      // Restore the state of the initial particles
      restoreParticles();

      // Delete newly created particles
      for( ParticleIter i = created.begin(); i != created.end(); ++i )
        delete *i;

      FinalState *fsBlocked = new FinalState;
      delete fs;
      fsBlocked->makeNoEnergyConservation();
      fsBlocked->setTotalEnergyBeforeInteraction(0.0);

      return fsBlocked; // Interaction is blocked. Return an empty final state.
    }
    DEBUG("Enforcing energy conservation: success!" << std::endl);

    // Check that outgoing delta resonances can decay to pi-N
    for( ParticleIter i = modified.begin(); i != modified.end(); ++i )
      if((*i)->isDelta() &&
          (*i)->getMass() < ParticleTable::effectiveDeltaDecayThreshold) {
        DEBUG("Mass of the produced delta below decay threshold; forbidding collision. deltaMass=" <<
            (*i)->getMass() << std::endl);

        // Restore the state of the initial particles
        restoreParticles();

        // Delete newly created particles
        for( ParticleIter j = created.begin(); j != created.end(); ++j )
          delete *j;

        FinalState *fsBlocked = new FinalState;
        delete fs;
        fsBlocked->makeNoEnergyConservation();
        fsBlocked->setTotalEnergyBeforeInteraction(0.0);

        return fsBlocked; // Interaction is blocked. Return an empty final state.
      }

    // Test Pauli blocking
    G4bool isBlocked = Pauli::isBlocked(modifiedAndCreated, theNucleus);

    if(isBlocked) {
      DEBUG("Pauli: Blocked!" << std::endl);

      // Restore the state of the initial particles
      restoreParticles();

      // Delete newly created particles
      for( ParticleIter i = created.begin(); i != created.end(); ++i )
        delete *i;

      FinalState *fsBlocked = new FinalState;
      delete fs;
      fsBlocked->makePauliBlocked();
      fsBlocked->setTotalEnergyBeforeInteraction(0.0);

      return fsBlocked; // Interaction is blocked. Return an empty final state.
    }
    DEBUG("Pauli: Allowed!" << std::endl);

    // Test CDPP blocking
    G4bool isCDPPBlocked = Pauli::isCDPPBlocked(created, theNucleus);

    if(isCDPPBlocked) {
      DEBUG("CDPP: Blocked!" << std::endl);

      // Restore the state of the initial particles
      restoreParticles();

      // Delete newly created particles
      for( ParticleIter i = created.begin(); i != created.end(); ++i )
        delete *i;

      FinalState *fsBlocked = new FinalState;
      delete fs;
      fsBlocked->makePauliBlocked();
      fsBlocked->setTotalEnergyBeforeInteraction(0.0);

      return fsBlocked; // Interaction is blocked. Return an empty final state.
    }
    DEBUG("CDPP: Allowed!" << std::endl);

    // If all went well, try to bring particles inside the nucleus...
    for( ParticleIter i = modifiedAndCreated.begin(); i != modifiedAndCreated.end(); ++i )
    {
      // ...except for pions beyond their surface radius.
      if((*i)->isOutOfWell()) continue;

      const G4bool successBringParticlesInside = bringParticleInside(*i);
      if( !successBringParticlesInside ) {
        ERROR("Failed to bring particle inside the nucleus!" << std::endl);
      }
    }

    // Collision accepted!
    for( ParticleIter i = modifiedAndCreated.begin(); i != modifiedAndCreated.end(); ++i ) {
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
          DEBUG("The following particle goes back to spectator:" << std::endl
              << (*i)->print() << std::endl);
          if(!(*i)->isTargetSpectator()) {
            theNucleus->getStore()->getBook()->decrementCascading();
          }
          (*i)->makeTargetSpectator();
        } else {
          if((*i)->isTargetSpectator()) {
            theNucleus->getStore()->getBook()->incrementCascading();
          }
          (*i)->makeParticipant();
        }
      }
    }
    ParticleList destroyed = fs->getDestroyedParticles();
    for( ParticleIter i = destroyed.begin(); i != destroyed.end(); ++i )
      if(!(*i)->isTargetSpectator())
        theNucleus->getStore()->getBook()->decrementCascading();

    return fs;
  }

  void InteractionAvatar::restoreParticles() const {
    particle1->setType(oldParticle1Type);
    particle1->setEnergy(oldParticle1Energy);
    particle1->setPotentialEnergy(oldParticle1Potential);
    particle1->setMomentum(oldParticle1Momentum);
    particle1->setPosition(oldParticle1Position);
    particle1->setMass(oldParticle1Mass);
    particle1->setHelicity(oldParticle1Helicity);

    if(particle2) {
      particle2->setType(oldParticle2Type);
      particle2->setEnergy(oldParticle2Energy);
      particle2->setPotentialEnergy(oldParticle2Potential);
      particle2->setMomentum(oldParticle2Momentum);
      particle2->setPosition(oldParticle2Position);
      particle2->setMass(oldParticle2Mass);
      particle2->setHelicity(oldParticle2Helicity);
    }
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
      violationEFunctor = new ViolationEEnergyFunctor(theNucleus, fs);
    }

    // Apply the root-finding algorithm
    const G4bool success = RootFinder::solve(violationEFunctor, 1.0);
    if(success) { // Apply the solution
      std::pair<G4double,G4double> theSolution = RootFinder::getSolution();
      (*violationEFunctor)(theSolution.first);
    } else {
      WARN("Couldn't enforce energy conservation after an interaction, root-finding algorithm failed." << std::endl);
    }
    delete violationEFunctor;
    return success;
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
    ParticleList created = finalState->getCreatedParticles();
    finalParticles.splice(finalParticles.end(), created);

    // Store the particle momenta (necessary for the calls to
    // scaleParticleMomenta() to work)
    particleMomenta.clear();
    for(ParticleIter i=finalParticles.begin(); i!=finalParticles.end(); ++i) {
      (*i)->boost(*boostVector);
      particleMomenta.push_back((*i)->getMomentum());
    }
  }

  G4double InteractionAvatar::ViolationEMomentumFunctor::operator()(const G4double alpha) const {
    scaleParticleMomenta(alpha);

    G4double deltaE = 0.0;
    for(ParticleIter i=finalParticles.begin(); i!=finalParticles.end(); ++i)
      deltaE += (*i)->getEnergy() - (*i)->getPotentialEnergy();
    deltaE -= initialEnergy;
    return deltaE;
  }

  void InteractionAvatar::ViolationEMomentumFunctor::scaleParticleMomenta(const G4double alpha) const {

    std::list<ThreeVector>::const_iterator iP = particleMomenta.begin();
    for(ParticleIter i=finalParticles.begin(); i!=finalParticles.end(); ++i, ++iP) {
      (*i)->setMomentum((*iP)*alpha);
      (*i)->adjustEnergyFromMomentum();
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

  InteractionAvatar::ViolationEEnergyFunctor::ViolationEEnergyFunctor(Nucleus * const nucleus, FinalState const * const finalState) :
    RootFunctor(0., 1E6),
    initialEnergy(finalState->getTotalEnergyBeforeInteraction()),
    theNucleus(nucleus),
    theParticle(finalState->getModifiedParticles().front()),
    theEnergy(theParticle->getEnergy()),
    theMomentum(theParticle->getMomentum()),
    energyThreshold(KinematicsUtils::energy(theMomentum,ParticleTable::effectiveDeltaDecayThreshold))
  {
// assert(theNucleus);
// assert(finalState->getModifiedParticles().size()==1);
// assert(theParticle->isDelta());
  }

  G4double InteractionAvatar::ViolationEEnergyFunctor::operator()(const G4double alpha) const {
    setParticleEnergy(alpha);
    return theParticle->getEnergy() - theParticle->getPotentialEnergy() - initialEnergy;
  }

  void InteractionAvatar::ViolationEEnergyFunctor::setParticleEnergy(const G4double alpha) const {

    G4double locE = KinematicsUtils::getLocalEnergy(theNucleus, theParticle); // Initial value of local energy
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
      theNucleus->updatePotentialEnergy(theParticle); // ...update its potential energy...
      locE = KinematicsUtils::getLocalEnergy(theNucleus, theParticle); // ...and recompute locE.
      deltaLocE = std::abs(locE-locEOld);
    }

  }

  void InteractionAvatar::ViolationEEnergyFunctor::cleanUp(const G4bool success) const {
    if(!success)
      setParticleEnergy(1.);
  }

}
