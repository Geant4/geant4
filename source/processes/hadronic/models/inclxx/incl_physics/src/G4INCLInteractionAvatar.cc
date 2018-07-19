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
    weight(1.),
    violationEFunctor(NULL)
  {
  }

  InteractionAvatar::InteractionAvatar(G4double time, G4INCL::Nucleus *n, G4INCL::Particle *p1,
      G4INCL::Particle *p2)
    : IAvatar(time), theNucleus(n),
    particle1(p1), particle2(p2),
    isPiN((p1->isPion() && p2->isNucleon()) || (p2->isPion() && p1->isNucleon())),
    weight(1.),
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
    if(!theNucleus || p->isMeson()) return; // Local energy does not make any sense without a nucleus

    if(shouldUseLocalEnergy())
      KinematicsUtils::transformToLocalEnergyFrame(theNucleus, p);
  }

  void InteractionAvatar::preInteraction() {
    preInteractionBlocking();

    preInteractionLocalEnergy(particle1);

    if(particle2) {
      preInteractionLocalEnergy(particle2);
      boostVector = KinematicsUtils::makeBoostVector(particle1, particle2);
      particle2->boost(boostVector);
    } else {
      boostVector = particle1->getMomentum()/particle1->getEnergy();
    }
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

    while( pos2 >= r*r && iterations<maxIterations ) /* Loop checking, 10.07.2015, D.Mancusi */
    {
      pos *= std::sqrt(r*r*0.9801/pos2); // 0.9801 == 0.99*0.99
      pos2 = pos.mag2();
      iterations++;
    }
    if( iterations < maxIterations)
    {
      INCL_DEBUG("Particle position vector length was : " << p->getPosition().mag() << ", rescaled to: " << pos.mag() << '\n');
      p->setPosition(pos);
      return true;
    }
    else
      return false;
  }

  void InteractionAvatar::postInteraction(FinalState *fs) {
    INCL_DEBUG("postInteraction: final state: " << '\n' << fs->print() << '\n');
    modified = fs->getModifiedParticles();
    created = fs->getCreatedParticles();
    Destroyed = fs->getDestroyedParticles();
    modifiedAndCreated = modified;
    modifiedAndCreated.insert(modifiedAndCreated.end(), created.begin(), created.end());
    ModifiedAndDestroyed = modified;
    ModifiedAndDestroyed.insert(ModifiedAndDestroyed.end(), Destroyed.begin(), Destroyed.end());

    // Boost back to lab
    modifiedAndCreated.boost(-boostVector);

    // If there is no Nucleus, just return
    if(!theNucleus) return;

    // Mark pions that have been created outside their well (we will force them
    // to be emitted later).
    for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i )
      if((*i)->isPion() && (*i)->getPosition().mag() > theNucleus->getSurfaceRadius(*i)) {
        (*i)->makeParticipant();
        (*i)->setOutOfWell();
        fs->addOutgoingParticle(*i);
        INCL_DEBUG("Pion was created outside its potential well." << '\n'
            << (*i)->print());
      }

    // Try to enforce energy conservation
    fs->setTotalEnergyBeforeInteraction(oldTotalEnergy);
    G4bool success = enforceEnergyConservation(fs);
    if(!success) {
      INCL_DEBUG("Enforcing energy conservation: failed!" << '\n');

      // Restore the state of the initial particles
      restoreParticles();

      // Delete newly created particles
      for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i )
        delete *i;

      fs->reset();
      fs->makeNoEnergyConservation();
      fs->setTotalEnergyBeforeInteraction(0.0);

      return; // Interaction is blocked. Return an empty final state.
    }
    INCL_DEBUG("Enforcing energy conservation: success!" << '\n');

    INCL_DEBUG("postInteraction after energy conservation: final state: " << '\n' << fs->print() << '\n');

    // Check that outgoing delta resonances can decay to pi-N
    for(ParticleIter i=modified.begin(), e=modified.end(); i!=e; ++i )
      if((*i)->isDelta() &&
          (*i)->getMass() < ParticleTable::minDeltaMass) {
        INCL_DEBUG("Mass of the produced delta below decay threshold; forbidding collision. deltaMass=" <<
            (*i)->getMass() << '\n');

        // Restore the state of the initial particles
        restoreParticles();

        // Delete newly created particles
        for(ParticleIter j=created.begin(), end=created.end(); j!=end; ++j )
          delete *j;

        fs->reset();
        fs->makeNoEnergyConservation();
        fs->setTotalEnergyBeforeInteraction(0.0);

        return; // Interaction is blocked. Return an empty final state.
      }

    INCL_DEBUG("Random seeds before Pauli blocking: " << Random::getSeeds() << '\n');
    // Test Pauli blocking
    G4bool isBlocked = Pauli::isBlocked(modifiedAndCreated, theNucleus);

    if(isBlocked) {
      INCL_DEBUG("Pauli: Blocked!" << '\n');

      // Restore the state of the initial particles
      restoreParticles();

      // Delete newly created particles
      for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i )
        delete *i;

      fs->reset();
      fs->makePauliBlocked();
      fs->setTotalEnergyBeforeInteraction(0.0);

      return; // Interaction is blocked. Return an empty final state.
    }
    INCL_DEBUG("Pauli: Allowed!" << '\n');

    // Test CDPP blocking
    G4bool isCDPPBlocked = Pauli::isCDPPBlocked(created, theNucleus);

    if(isCDPPBlocked) {
      INCL_DEBUG("CDPP: Blocked!" << '\n');

      // Restore the state of the initial particles
      restoreParticles();

      // Delete newly created particles
      for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i )
        delete *i;

      fs->reset();
      fs->makePauliBlocked();
      fs->setTotalEnergyBeforeInteraction(0.0);

      return; // Interaction is blocked. Return an empty final state.
    }
    INCL_DEBUG("CDPP: Allowed!" << '\n');

    // If all went well, try to bring particles inside the nucleus...
    for(ParticleIter i=modifiedAndCreated.begin(), e=modifiedAndCreated.end(); i!=e; ++i )
    {
      // ...except for pions beyond their surface radius.
      if((*i)->isOutOfWell()) continue;

      const G4bool successBringParticlesInside = bringParticleInside(*i);
      if( !successBringParticlesInside ) {
        INCL_ERROR("Failed to bring particle inside the nucleus!" << '\n');
      }
    }

    // Collision accepted!
    // Biasing of the final state
    std::vector<G4int> newBiasCollisionVector;
    newBiasCollisionVector = ModifiedAndDestroyed.getParticleListBiasVector();
    if(std::fabs(weight-1.) > 1E-6){
	  newBiasCollisionVector.push_back(Particle::nextBiasedCollisionID);
	  Particle::FillINCLBiasVector(1./weight);
	  weight = 1.; //Should be reinitialized in case of next collision non baised
	}
    for(ParticleIter i=modifiedAndCreated.begin(), e=modifiedAndCreated.end(); i!=e; ++i ) {
	  (*i)->setBiasCollisionVector(newBiasCollisionVector);
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
          INCL_DEBUG("The following particle goes back to spectator:" << '\n'
              << (*i)->print() << '\n');
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
    return;
  }

  void InteractionAvatar::restoreParticles() const {
    (*particle1) = (*backupParticle1);
    if(particle2)
      (*particle2) = (*backupParticle2);
  }

  G4bool InteractionAvatar::shouldUseLocalEnergy() const {
    if(!theNucleus) return false;
    LocalEnergyType theLocalEnergyType;
    if(getType()==DecayAvatarType || isPiN)
      theLocalEnergyType = theNucleus->getStore()->getConfig()->getLocalEnergyPiType();
    else
      theLocalEnergyType = theNucleus->getStore()->getConfig()->getLocalEnergyBBType();

    const G4bool firstAvatar = (theNucleus->getStore()->getBook().getAcceptedCollisions() == 0);
    return ((theLocalEnergyType == FirstCollisionLocalEnergy && firstAvatar) ||
            theLocalEnergyType == AlwaysLocalEnergy);
  }

  G4bool InteractionAvatar::enforceEnergyConservation(FinalState * const fs) {
    // Set up the violationE calculation
    const G4bool manyBodyFinalState = (modifiedAndCreated.size() > 1);
    
    if(manyBodyFinalState)
      violationEFunctor = new ViolationEMomentumFunctor(theNucleus, modifiedAndCreated, fs->getTotalEnergyBeforeInteraction(), boostVector, shouldUseLocalEnergy());
    else {
      Particle * const p = modified.front();
      // The following condition is necessary for the functor to work
      // correctly. A similar condition exists in INCL4.6.
      if(p->getMass() < ParticleTable::minDeltaMass)
        return false;
      violationEFunctor = new ViolationEEnergyFunctor(theNucleus, p, fs->getTotalEnergyBeforeInteraction(), shouldUseLocalEnergy());
    }

    // Apply the root-finding algorithm
    const RootFinder::Solution theSolution = RootFinder::solve(violationEFunctor, 1.0);
    if(theSolution.success) { // Apply the solution
      (*violationEFunctor)(theSolution.x);
    } else if(theNucleus){
      INCL_DEBUG("Couldn't enforce energy conservation after an interaction, root-finding algorithm failed." << '\n');
      theNucleus->getStore()->getBook().incrementEnergyViolationInteraction();
    }
    delete violationEFunctor;
    violationEFunctor = NULL;
    return theSolution.success;
  }

  /* ***                                                      ***
   * *** InteractionAvatar::ViolationEMomentumFunctor methods ***
   * ***                                                      ***/

  InteractionAvatar::ViolationEMomentumFunctor::ViolationEMomentumFunctor(Nucleus * const nucleus, ParticleList const &modAndCre, const G4double totalEnergyBeforeInteraction, ThreeVector const &boost, const G4bool localE) :
    RootFunctor(0., 1E6),
    finalParticles(modAndCre),
    initialEnergy(totalEnergyBeforeInteraction),
    theNucleus(nucleus),
    boostVector(boost),
    shouldUseLocalEnergy(localE)
  {
    // Store the particle momenta (necessary for the calls to
    // scaleParticleMomenta() to work)
    for(ParticleIter i=finalParticles.begin(), e=finalParticles.end(); i!=e; ++i) {
      (*i)->boost(boostVector);
      particleMomenta.push_back((*i)->getMomentum());
    }
  }

  InteractionAvatar::ViolationEMomentumFunctor::~ViolationEMomentumFunctor() {
    particleMomenta.clear();
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

    std::vector<ThreeVector>::const_iterator iP = particleMomenta.begin();
    for(ParticleIter i=finalParticles.begin(), e=finalParticles.end(); i!=e; ++i, ++iP) {
      (*i)->setMomentum((*iP)*alpha);
      (*i)->adjustEnergyFromMomentum();
      (*i)->rpCorrelate();
      (*i)->boost(-boostVector);
      if(theNucleus)
        theNucleus->updatePotentialEnergy(*i);
      else
        (*i)->setPotentialEnergy(0.);

//jcd      if(shouldUseLocalEnergy && !(*i)->isPion()) { // This translates AECSVT's loops 1, 3 and 4
            if(shouldUseLocalEnergy && !(*i)->isPion() && !(*i)->isEta() && !(*i)->isOmega() &&
                                       !(*i)->isKaon() && !(*i)->isAntiKaon() &&!(*i)->isLambda() && !(*i)->isSigma()) { // This translates AECSVT's loops 1, 3 and 4
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

  InteractionAvatar::ViolationEEnergyFunctor::ViolationEEnergyFunctor(Nucleus * const nucleus, Particle * const aParticle, const G4double totalEnergyBeforeInteraction, const G4bool localE) :
    RootFunctor(0., 1E6),
    initialEnergy(totalEnergyBeforeInteraction),
    theNucleus(nucleus),
    theParticle(aParticle),
    theEnergy(theParticle->getEnergy()),
    theMomentum(theParticle->getMomentum()),
    energyThreshold(KinematicsUtils::energy(theMomentum,ParticleTable::minDeltaMass)),
    shouldUseLocalEnergy(localE)
  {
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
      G4double particleEnergy = energyThreshold + locE + alpha*(theEnergy-energyThreshold);
      const G4double theMass2 = std::pow(particleEnergy,2.)-theMomentum.mag2();
      G4double theMass;
      if(theMass2>ParticleTable::minDeltaMass2)
        theMass = std::sqrt(theMass2);
      else {
        theMass = ParticleTable::minDeltaMass;
        particleEnergy = energyThreshold;
      }
      theParticle->setMass(theMass);
      theParticle->setEnergy(particleEnergy); // Update the energy of the particle...
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
