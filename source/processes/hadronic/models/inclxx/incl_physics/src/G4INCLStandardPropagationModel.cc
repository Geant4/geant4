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

/*
 * StandardPropagationModel.cpp
 *
 *  \date 4 juin 2009
 * \author Pekka Kaitaniemi
 */

#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLSurfaceAvatar.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLDecayAvatar.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLRandom.hh"
#include <iostream>
#include <algorithm>
#include "G4INCLLogger.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLCoulombDistortion.hh"
#include "G4INCLDeltaDecayChannel.hh"
#include "G4INCLParticleEntryAvatar.hh"
#include "G4INCLIntersection.hh"

namespace G4INCL {

    StandardPropagationModel::StandardPropagationModel(LocalEnergyType localEnergyType, LocalEnergyType localEnergyDeltaType)
      :theNucleus(0), maximumTime(70.0), currentTime(0.0), firstAvatar(true),
      theLocalEnergyType(localEnergyType),
      theLocalEnergyDeltaType(localEnergyDeltaType)
    {
    }

    StandardPropagationModel::~StandardPropagationModel()
    {
      delete theNucleus;
    }

    G4INCL::Nucleus* StandardPropagationModel::getNucleus()
    {
      return theNucleus;
    }

    G4double StandardPropagationModel::shoot(ParticleSpecies const projectileSpecies, const G4double kineticEnergy, const G4double impactParameter, const G4double phi) {
      if(projectileSpecies.theType==Composite)
        return shootComposite(projectileSpecies, kineticEnergy, impactParameter, phi);
      else
        return shootParticle(projectileSpecies.theType, kineticEnergy, impactParameter, phi);
    }

    G4double StandardPropagationModel::shootParticle(ParticleType const type, const G4double kineticEnergy, const G4double impactParameter, const G4double phi) {
      theNucleus->setParticleNucleusCollision();
      currentTime = 0.0;

      // Create the projectile particle
      const G4double projectileMass = ParticleTable::getTableParticleMass(type);
      G4double energy = kineticEnergy + projectileMass;
      G4double momentumZ = std::sqrt(energy*energy - projectileMass*projectileMass);
      ThreeVector momentum(0.0, 0.0, momentumZ);
      Particle *p= new G4INCL::Particle(type, energy, momentum, ThreeVector());

      G4double temfin = 0.0;
      if( p->isNucleon() )
        temfin = 29.8 * std::pow(theNucleus->getA(), 0.16);
      else {
        const G4double tlab = p->getEnergy() - p->getMass();
        temfin = 30.18 * std::pow(theNucleus->getA(), 0.17*(1.0 - 5.7E-5*tlab));
      }

      maximumTime = temfin;

      // If the incoming particle is slow, use a larger stopping time
      const G4double rMax = theNucleus->getUniverseRadius();
      const G4double distance = 2.*rMax;
      const G4double projectileVelocity = p->boostVector().mag();
      const G4double traversalTime = distance / projectileVelocity;
      if(maximumTime < traversalTime)
        maximumTime = traversalTime;
      DEBUG("Cascade stopping time is " << maximumTime << std::endl);

      // If Coulomb is activated, do not process events with impact
      // parameter larger than the maximum impact parameter, taking into
      // account Coulomb distortion.
      if(impactParameter>CoulombDistortion::maxImpactParameter(p->getSpecies(), kineticEnergy, theNucleus)) {
        DEBUG("impactParameter>CoulombDistortion::maxImpactParameter" << std::endl);
        delete p;
        return -1.;
      }

      ThreeVector position(impactParameter * std::cos(phi),
          impactParameter * std::sin(phi),
          0.);
      p->setPosition(position);

      // Fill in the relevant kinematic variables
      theNucleus->setIncomingAngularMomentum(p->getAngularMomentum());
      theNucleus->setIncomingMomentum(p->getMomentum());
      theNucleus->setInitialEnergy(p->getEnergy()
          + ParticleTable::getTableMass(theNucleus->getA(),theNucleus->getZ()));

      // Reset the particle kinematics to the INCL values
      p->setINCLMass();
      p->setEnergy(p->getMass() + kineticEnergy);
      p->adjustMomentumFromEnergy();

      p->makeProjectileSpectator();
      generateAllAvatars();
      firstAvatar = false;

      // Get the entry avatars from Coulomb and put them in the Store
      ParticleEntryAvatar *theEntryAvatar = CoulombDistortion::bringToSurface(p, theNucleus);
      if(theEntryAvatar) {
        theNucleus->getStore()->addParticleEntryAvatar(theEntryAvatar);

        theNucleus->setProjectileChargeNumber(p->getZ());
        theNucleus->setProjectileMassNumber(p->getA());

        return p->getTransversePosition().mag();
      } else {
        delete p;
        return -1.;
      }
    }

    G4double StandardPropagationModel::shootComposite(ParticleSpecies const species, const G4double kineticEnergy, const G4double impactParameter, const G4double phi) {
      theNucleus->setNucleusNucleusCollision();
      currentTime = 0.0;

      // Create the ProjectileRemnant object
      ProjectileRemnant *pr = new ProjectileRemnant(species, kineticEnergy);

      // Same stopping time as for nucleon-nucleus
      maximumTime = 29.8 * std::pow(theNucleus->getA(), 0.16);

      // If the incoming cluster is slow, use a larger stopping time
      const G4double rms = ParticleTable::getNuclearRadius(pr->getA(), pr->getZ());
      const G4double rMax = theNucleus->getUniverseRadius();
      const G4double distance = 2.*rMax + 2.725*rms;
      const G4double projectileVelocity = pr->boostVector().mag();
      const G4double traversalTime = distance / projectileVelocity;
      if(maximumTime < traversalTime)
        maximumTime = traversalTime;
      DEBUG("Cascade stopping time is " << maximumTime << std::endl);

      // If Coulomb is activated, do not process events with impact
      // parameter larger than the maximum impact parameter, taking into
      // account Coulomb distortion.
      if(impactParameter>CoulombDistortion::maxImpactParameter(pr,theNucleus)) {
        pr->deleteParticles();
        DEBUG("impactParameter>CoulombDistortion::maxImpactParameter" << std::endl);
        delete pr;
        return -1.;
      }

      // Position the cluster at the right impact parameter
      ThreeVector position(impactParameter * std::cos(phi),
          impactParameter * std::sin(phi),
          0.);
      pr->setPosition(position);

      /* Store the internal kinematics of the projectile remnant.
       *
       * Note that this is at variance with the Fortran version, which stores
       * the initial kinematics of the particles *after* putting the spectators
       * on mass shell, but *before* removing the same energy from the
       * participants. Due to the different code flow, doing so in the C++
       * version leads to wrong excitation energies for the forced compound
       * nucleus.
       */
      pr->storeComponents();

      // Fill in the relevant kinematic variables
      theNucleus->setIncomingAngularMomentum(pr->getAngularMomentum());
      theNucleus->setIncomingMomentum(pr->getMomentum());
      theNucleus->setInitialEnergy(pr->getEnergy()
          + ParticleTable::getTableMass(theNucleus->getA(),theNucleus->getZ()));

      generateAllAvatars();
      firstAvatar = false;

      // Get the entry avatars from Coulomb
      IAvatarList theAvatarList
        = CoulombDistortion::bringToSurface(pr, theNucleus);

      if(theAvatarList.empty()) {
        DEBUG("No ParticleEntryAvatar found, transparent event" << std::endl);
        pr->deleteParticles();
        delete pr;
        return -1.;
      }

      // Tell the Nucleus about the ProjectileRemnant
      theNucleus->setProjectileRemnant(pr);

      // Set the number of projectile particles
      theNucleus->setProjectileChargeNumber(pr->getZ());
      theNucleus->setProjectileMassNumber(pr->getA());

      // Register the ParticleEntryAvatars
      theNucleus->getStore()->addParticleEntryAvatars(theAvatarList);

      return pr->getTransversePosition().mag();
    }

    G4double StandardPropagationModel::getStoppingTime() {
      return maximumTime;
    }

    void StandardPropagationModel::setStoppingTime(G4double time) {
// assert(time>0.0);
      maximumTime = time;
    }

    G4double StandardPropagationModel::getCurrentTime() {
      return currentTime;
    }

    void StandardPropagationModel::setNucleus(G4INCL::Nucleus *nucleus)
    {
      theNucleus = nucleus;
    }

    void StandardPropagationModel::registerAvatar(G4INCL::IAvatar *anAvatar)
    {
      if(anAvatar) theNucleus->getStore()->add(anAvatar);
    }

    IAvatar *StandardPropagationModel::generateBinaryCollisionAvatar(Particle * const p1, Particle * const p2) const {
      // Is either particle a participant?
      if(!p1->isParticipant() && !p2->isParticipant() && p1->getParticipantType()==p2->getParticipantType()) return NULL;

      // Is it a pi-resonance collision (we don't treat them)?
      if((p1->isResonance() && p2->isPion()) || (p1->isPion() && p2->isResonance()))
        return NULL;

      // Will the avatar take place between now and the end of the cascade?
      G4double minDistOfApproachSquared = 0.0;
      G4double t = getTime(p1, p2, &minDistOfApproachSquared);
      if(t>maximumTime || t<currentTime) return NULL;

      // Local energy. Jump through some hoops to calculate the cross section
      // at the collision point, and clean up after yourself afterwards.
      G4bool hasLocalEnergy;
      if(p1->isPion() || p2->isPion())
        hasLocalEnergy = ((theLocalEnergyDeltaType == FirstCollisionLocalEnergy &&
              theNucleus->getStore()->getBook()->getAcceptedCollisions()==0) ||
            theLocalEnergyDeltaType == AlwaysLocalEnergy);
      else
        hasLocalEnergy = ((theLocalEnergyType == FirstCollisionLocalEnergy &&
              theNucleus->getStore()->getBook()->getAcceptedCollisions()==0) ||
            theLocalEnergyType == AlwaysLocalEnergy);
      const G4bool p1HasLocalEnergy = (hasLocalEnergy && !p1->isPion());
      const G4bool p2HasLocalEnergy = (hasLocalEnergy && !p2->isPion());

      Particle backupParticle1 = *p1;
      if(p1HasLocalEnergy) {
        p1->propagate(t - currentTime);
        if(p1->getPosition().mag() > theNucleus->getSurfaceRadius(p1)) {
          *p1 = backupParticle1;
          return NULL;
        }
        KinematicsUtils::transformToLocalEnergyFrame(theNucleus, p1);
      }
      Particle backupParticle2 = *p2;
      if(p2HasLocalEnergy) {
        p2->propagate(t - currentTime);
        if(p2->getPosition().mag() > theNucleus->getSurfaceRadius(p2)) {
          *p2 = backupParticle2;
          if(p1HasLocalEnergy) {
            *p1 = backupParticle1;
          }
          return NULL;
        }
        KinematicsUtils::transformToLocalEnergyFrame(theNucleus, p2);
      }

      // Compute the total cross section
      const G4double totalCrossSection = CrossSections::total(p1, p2);
      const G4double squareTotalEnergyInCM = KinematicsUtils::squareTotalEnergyInCM(p1,p2);

      // Restore particles to their state before the local-energy tweak
      if(p1HasLocalEnergy) {
        *p1 = backupParticle1;
      }
      if(p2HasLocalEnergy) {
        *p2 = backupParticle2;
      }

      // Is the CM energy > cutNN? (no cutNN on the first collision)
      if(theNucleus->getStore()->getBook()->getAcceptedCollisions()>0
          && p1->isNucleon() && p2->isNucleon()
          && squareTotalEnergyInCM < BinaryCollisionAvatar::cutNNSquared) return NULL;

      // Do the particles come close enough to each other?
      if(Math::tenPi*minDistOfApproachSquared > totalCrossSection) return NULL;

      // Bomb out if the two collision partners are the same particle
// assert(p1->getID() != p2->getID());

      // Return a new avatar, then!
      return new G4INCL::BinaryCollisionAvatar(t, totalCrossSection, theNucleus, p1, p2);
    }

    G4double StandardPropagationModel::getReflectionTime(G4INCL::Particle const * const aParticle) {
      Intersection theIntersection(
          IntersectionFactory::getLaterTrajectoryIntersection(
            aParticle->getPosition(),
            aParticle->getPropagationVelocity(),
            theNucleus->getSurfaceRadius(aParticle)));
      G4double time;
      if(theIntersection.exists) {
        time = currentTime + theIntersection.time;
      } else {
        ERROR("Imaginary reflection time for particle: " << std::endl
          << aParticle->print());
        time = 10000.0;
      }
      return time;
    }

    G4double StandardPropagationModel::getTime(G4INCL::Particle const * const particleA,
					     G4INCL::Particle const * const particleB, G4double *minDistOfApproach) const
    {
      G4double time;
      G4INCL::ThreeVector t13 = particleA->getPropagationVelocity();
      t13 -= particleB->getPropagationVelocity();
      G4INCL::ThreeVector distance = particleA->getPosition();
      distance -= particleB->getPosition();
      const G4double t7 = t13.dot(distance);
      const G4double dt = t13.mag2();
      if(dt <= 1.0e-10) {
        (*minDistOfApproach) = 100000.0;
        return currentTime + 100000.0;
      } else {
        time = -t7/dt;
      }
      (*minDistOfApproach) = distance.mag2() + time * t7;
      return currentTime + time;
    }

    void StandardPropagationModel::generateUpdatedCollisions(const ParticleList &updatedParticles, const ParticleList &particles) {

      // Loop over all the updated particles
      for(ParticleIter updated = updatedParticles.begin(); updated != updatedParticles.end(); ++updated)
      {
        // Loop over all the particles
        for(ParticleIter particle = particles.begin(); particle != particles.end(); ++particle)
        {
          /* Consider the generation of a collision avatar only if (*particle)
           * is not one of the updated particles.
           * The criterion makes sure that you don't generate avatars between
           * updated particles. */
          if((*particle)->isInList(updatedParticles)) continue;

          registerAvatar(generateBinaryCollisionAvatar(*particle,*updated));
        }
      }
    }

    void StandardPropagationModel::generateCollisions(const ParticleList &particles, const ParticleList &except) {

      G4bool haveExcept;
      haveExcept=(except.size()!=0);

      // Loop over all the particles
      for(ParticleIter p1 = particles.begin(); p1 != particles.end(); ++p1)
      {
        // Loop over the rest of the particles
        ParticleIter p2 = p1;
        for(++p2; p2 != particles.end(); ++p2)
        {
          // Skip the collision if both particles must be excluded
          if(haveExcept && (*p1)->isInList(except) && (*p2)->isInList(except)) continue;

          registerAvatar(generateBinaryCollisionAvatar(*p1,*p2));
        }
      }

    }

    void StandardPropagationModel::updateAvatars(const ParticleList &particles) {

      for(ParticleIter iter = particles.begin(); iter != particles.end(); ++iter) {
        G4double time = this->getReflectionTime(*iter);
        if(time <= maximumTime) registerAvatar(new SurfaceAvatar(*iter, time, theNucleus));
      }
      ParticleList const &p = theNucleus->getStore()->getParticles();
      generateUpdatedCollisions(particles, p);             // Predict collisions with spectators and participants
    }

    void StandardPropagationModel::generateAllAvatars(G4bool excludeUpdated) {
      ParticleList particles = theNucleus->getStore()->getParticles();
      if(particles.empty()) { ERROR("No particles inside the nucleus!" << std::endl); }
      for(ParticleIter i = particles.begin(); i != particles.end(); ++i) {
        G4double time = this->getReflectionTime(*i);
        if(time <= maximumTime) registerAvatar(new SurfaceAvatar(*i, time, theNucleus));
      }
      ParticleList except;
      if(excludeUpdated)
        except = theNucleus->getUpdatedParticles();
      generateCollisions(particles,except);
      generateDecays(particles);
    }

    void StandardPropagationModel::generateDecays(const ParticleList &particles) {
      for(ParticleIter i = particles.begin(); i != particles.end(); ++i) {
	if((*i)->isDelta()) {
    G4double decayTime = DeltaDecayChannel::computeDecayTime((*i));
	  G4double time = currentTime + decayTime;
	  if(time <= maximumTime) {
	    registerAvatar(new DecayAvatar((*i), time, theNucleus));
	  }
	}
      }
    }

    G4INCL::IAvatar* StandardPropagationModel::propagate()
    {
      // We update only the information related to particles that were updated
      // by the previous avatar.
#ifdef INCL_REGENERATE_AVATARS
#warning "The INCL_REGENERATE_AVATARS code has not been tested in a while. Use it at your peril."
      if(theNucleus->getUpdatedParticles().size()!=0 || theNucleus->getCreatedParticles().size()!=0) {
        // Regenerates the entire avatar list, skipping collisions between
        // updated particles
        theNucleus->getStore()->clearAvatars();
        theNucleus->getStore()->initialiseParticleAvatarConnections();
        generateAllAvatars(true);
      }
#else
      // Deltas are created by transforming nucleon into a delta for
      // efficiency reasons
      Particle * const blockedDelta = theNucleus->getBlockedDelta();
      ParticleList updatedParticles = theNucleus->getUpdatedParticles();
      if(blockedDelta)
        updatedParticles.push_back(blockedDelta);
      generateDecays(updatedParticles);

      ParticleList needNewAvatars = theNucleus->getUpdatedParticles();
      ParticleList created = theNucleus->getCreatedParticles();
      needNewAvatars.splice(needNewAvatars.end(), created);
      updateAvatars(needNewAvatars);
#endif

      G4INCL::IAvatar *theAvatar = theNucleus->getStore()->findSmallestTime();
      if(theAvatar == 0) return 0; // Avatar list is empty
      //      theAvatar->dispose();

      if(theAvatar->getTime() < currentTime) {
        ERROR("Avatar time = " << theAvatar->getTime() << ", currentTime = " << currentTime << std::endl);
        return 0;
      } else if(theAvatar->getTime() > currentTime) {
        theNucleus->getStore()->timeStep(theAvatar->getTime() - currentTime);

        currentTime = theAvatar->getTime();
        theNucleus->getStore()->getBook()->setCurrentTime(currentTime);
      }

      return theAvatar;
    }

    void StandardPropagationModel::putSpectatorsOnShell(IAvatarList const &entryAvatars, ParticleList const &spectators) {
      G4double deltaE = 0.0;
      for(ParticleIter p=spectators.begin(); p!=spectators.end(); ++p) {
        // put the spectators on shell (conserving their momentum)
        const G4double oldEnergy = (*p)->getEnergy();
        (*p)->setTableMass();
        (*p)->adjustEnergyFromMomentum();
        deltaE += (*p)->getEnergy() - oldEnergy;
      }

      deltaE /= entryAvatars.size(); // energy to remove from each participant

      for(IAvatarIter a=entryAvatars.begin(); a!=entryAvatars.end(); ++a) {
        // remove the energy from the participant
        Particle *p = (*a)->getParticles().front();
        ParticleType const t = p->getType();
        // we also need to slightly correct the participant mass
        const G4double energy = p->getEnergy() - deltaE
          - ParticleTable::getTableParticleMass(t) + ParticleTable::getINCLMass(t);
        p->setEnergy(energy);
        const G4double newMass = std::sqrt(energy*energy - p->getMomentum().mag2());
        p->setMass(newMass);
      }
    }
}
