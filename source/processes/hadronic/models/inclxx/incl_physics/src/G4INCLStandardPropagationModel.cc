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
#include "G4INCLSigmaZeroDecayChannel.hh"
#include "G4INCLPionResonanceDecayChannel.hh"
#include "G4INCLParticleEntryAvatar.hh"
#include "G4INCLIntersection.hh"

namespace G4INCL {

    StandardPropagationModel::StandardPropagationModel(LocalEnergyType localEnergyType, LocalEnergyType localEnergyDeltaType, const G4double hTime)
      :theNucleus(0), maximumTime(70.0), currentTime(0.0),
      hadronizationTime(hTime),
      firstAvatar(true),
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

    G4double StandardPropagationModel::shoot(ParticleSpecies const &projectileSpecies, const G4double kineticEnergy, const G4double impactParameter, const G4double phi) {
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

      G4double temfin;
      G4double TLab;
      if( p->isMeson()) {
        temfin = 30.18 * std::pow(theNucleus->getA(), 0.17);
        TLab = p->getKineticEnergy();
      } else {
        temfin = 29.8 * std::pow(theNucleus->getA(), 0.16);
        TLab = p->getKineticEnergy()/p->getA();
      }

      // energy-dependent stopping time above 2 AGeV
      if(TLab>2000.)
        temfin *= (5.8E4-TLab)/5.6E4;

      maximumTime = temfin;

      // If the incoming particle is slow, use a larger stopping time
      const G4double rMax = theNucleus->getUniverseRadius();
      const G4double distance = 2.*rMax;
      const G4double projectileVelocity = p->boostVector().mag();
      const G4double traversalTime = distance / projectileVelocity;
      if(maximumTime < traversalTime)
        maximumTime = traversalTime;
      INCL_DEBUG("Cascade stopping time is " << maximumTime << '\n');

      // If Coulomb is activated, do not process events with impact
      // parameter larger than the maximum impact parameter, taking into
      // account Coulomb distortion.
      if(impactParameter>CoulombDistortion::maxImpactParameter(p->getSpecies(), kineticEnergy, theNucleus)) {
        INCL_DEBUG("impactParameter>CoulombDistortion::maxImpactParameter" << '\n');
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

        return p->getTransversePosition().mag();
      } else {
        delete p;
        return -1.;
      }
    }

    G4double StandardPropagationModel::shootComposite(ParticleSpecies const &species, const G4double kineticEnergy, const G4double impactParameter, const G4double phi) {
      theNucleus->setNucleusNucleusCollision();
      currentTime = 0.0;

      // Create the ProjectileRemnant object
      ProjectileRemnant *pr = new ProjectileRemnant(species, kineticEnergy);

      // Same stopping time as for nucleon-nucleus
      maximumTime = 29.8 * std::pow(theNucleus->getA(), 0.16);

      // If the incoming cluster is slow, use a larger stopping time
      const G4double rms = ParticleTable::getLargestNuclearRadius(pr->getA(), pr->getZ());
      const G4double rMax = theNucleus->getUniverseRadius();
      const G4double distance = 2.*rMax + 2.725*rms;
      const G4double projectileVelocity = pr->boostVector().mag();
      const G4double traversalTime = distance / projectileVelocity;
      if(maximumTime < traversalTime)
        maximumTime = traversalTime;
      INCL_DEBUG("Cascade stopping time is " << maximumTime << '\n');

      // If Coulomb is activated, do not process events with impact
      // parameter larger than the maximum impact parameter, taking into
      // account Coulomb distortion.
      if(impactParameter>CoulombDistortion::maxImpactParameter(pr,theNucleus)) {
        INCL_DEBUG("impactParameter>CoulombDistortion::maxImpactParameter" << '\n');
        delete pr;
        return -1.;
      }

      // Position the cluster at the right impact parameter
      ThreeVector position(impactParameter * std::cos(phi),
          impactParameter * std::sin(phi),
          0.);
      pr->setPosition(position);

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
        INCL_DEBUG("No ParticleEntryAvatar found, transparent event" << '\n');
        delete pr;
        return -1.;
      }

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

      // Tell the Nucleus about the ProjectileRemnant
      theNucleus->setProjectileRemnant(pr);

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

    IAvatar *StandardPropagationModel::generateBinaryCollisionAvatar(Particle * const p1, Particle * const p2) {
      // Is either particle a participant?
      if(!p1->isParticipant() && !p2->isParticipant() && p1->getParticipantType()==p2->getParticipantType()) return NULL;

      // Is it a pi-resonance collision (we don't treat them)?
      if((p1->isResonance() && p2->isPion()) || (p1->isPion() && p2->isResonance()))
        return NULL;

      // Will the avatar take place between now and the end of the cascade?
      G4double minDistOfApproachSquared = 0.0;
      G4double t = getTime(p1, p2, &minDistOfApproachSquared);
      if(t>maximumTime || t<currentTime+hadronizationTime) return NULL;

      // Local energy. Jump through some hoops to calculate the cross section
      // at the collision point, and clean up after yourself afterwards.
      G4bool hasLocalEnergy;
      if(p1->isPion() || p2->isPion())
        hasLocalEnergy = ((theLocalEnergyDeltaType == FirstCollisionLocalEnergy &&
              theNucleus->getStore()->getBook().getAcceptedCollisions()==0) ||
            theLocalEnergyDeltaType == AlwaysLocalEnergy);
      else
        hasLocalEnergy = ((theLocalEnergyType == FirstCollisionLocalEnergy &&
              theNucleus->getStore()->getBook().getAcceptedCollisions()==0) ||
            theLocalEnergyType == AlwaysLocalEnergy);
      const G4bool p1HasLocalEnergy = (hasLocalEnergy && !p1->isMeson());
      const G4bool p2HasLocalEnergy = (hasLocalEnergy && !p2->isMeson());

      if(p1HasLocalEnergy) {
        backupParticle1 = *p1;
        p1->propagate(t - currentTime);
        if(p1->getPosition().mag() > theNucleus->getSurfaceRadius(p1)) {
          *p1 = backupParticle1;
          return NULL;
        }
        KinematicsUtils::transformToLocalEnergyFrame(theNucleus, p1);
      }
      if(p2HasLocalEnergy) {
        backupParticle2 = *p2;
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
      if(theNucleus->getStore()->getBook().getAcceptedCollisions()>0
          && p1->isNucleon() && p2->isNucleon()
          && squareTotalEnergyInCM < BinaryCollisionAvatar::getCutNNSquared()) return NULL;

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
        INCL_ERROR("Imaginary reflection time for particle: " << '\n'
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
      for(ParticleIter updated=updatedParticles.begin(), e=updatedParticles.end(); updated!=e; ++updated)
      {
        // Loop over all the particles
        for(ParticleIter particle=particles.begin(), end=particles.end(); particle!=end; ++particle)
        {
          /* Consider the generation of a collision avatar only if (*particle)
           * is not one of the updated particles.
           * The criterion makes sure that you don't generate avatars between
           * updated particles. */
          if(updatedParticles.contains(*particle)) continue;

          registerAvatar(generateBinaryCollisionAvatar(*particle,*updated));
        }
      }
    }

    void StandardPropagationModel::generateCollisions(const ParticleList &particles) {
      // Loop over all the particles
      for(ParticleIter p1=particles.begin(), e=particles.end(); p1!=e; ++p1) {
        // Loop over the rest of the particles
        for(ParticleIter p2 = p1 + 1; p2 != particles.end(); ++p2) {
          registerAvatar(generateBinaryCollisionAvatar(*p1,*p2));
        }
      }
    }

    void StandardPropagationModel::generateCollisions(const ParticleList &particles, const ParticleList &except) {

      const G4bool haveExcept = !except.empty();

      // Loop over all the particles
      for(ParticleIter p1=particles.begin(), e=particles.end(); p1!=e; ++p1)
      {
        // Loop over the rest of the particles
        ParticleIter p2 = p1;
        for(++p2; p2 != particles.end(); ++p2)
        {
          // Skip the collision if both particles must be excluded
          if(haveExcept && except.contains(*p1) && except.contains(*p2)) continue;

          registerAvatar(generateBinaryCollisionAvatar(*p1,*p2));
        }
      }

    }

    void StandardPropagationModel::updateAvatars(const ParticleList &particles) {

      for(ParticleIter iter=particles.begin(), e=particles.end(); iter!=e; ++iter) {
        G4double time = this->getReflectionTime(*iter);
        if(time <= maximumTime) registerAvatar(new SurfaceAvatar(*iter, time, theNucleus));
      }
      ParticleList const &p = theNucleus->getStore()->getParticles();
      generateUpdatedCollisions(particles, p);             // Predict collisions with spectators and participants
    }

    void StandardPropagationModel::generateAllAvatars() {
      ParticleList const &particles = theNucleus->getStore()->getParticles();
// assert(!particles.empty());
      G4double time;
      for(ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
        time = this->getReflectionTime(*i);
        if(time <= maximumTime) registerAvatar(new SurfaceAvatar(*i, time, theNucleus));
      }
      generateCollisions(particles);
      generateDecays(particles);
    }

#ifdef INCL_REGENERATE_AVATARS
    void StandardPropagationModel::generateAllAvatarsExceptUpdated(FinalState const * const fs) {
      ParticleList const &particles = theNucleus->getStore()->getParticles();
// assert(!particles.empty());
      for(ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
        G4double time = this->getReflectionTime(*i);
        if(time <= maximumTime) registerAvatar(new SurfaceAvatar(*i, time, theNucleus));
      }
      ParticleList except = fs->getModifiedParticles();
      ParticleList const &entering = fs->getEnteringParticles();
      except.insert(except.end(), entering.begin(), entering.end());
      generateCollisions(particles,except);
      generateDecays(particles);
    }
#endif

    void StandardPropagationModel::generateDecays(const ParticleList &particles) {
      for(ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
	       if((*i)->isDelta()) {
             G4double decayTime = DeltaDecayChannel::computeDecayTime((*i)); // time in fm/c
	         G4double time = currentTime + decayTime;
	         if(time <= maximumTime) {
	           registerAvatar(new DecayAvatar((*i), time, theNucleus));
	         }
	       }
	       else if((*i)->getType() == SigmaZero) {
	         G4double decayTime = SigmaZeroDecayChannel::computeDecayTime((*i)); // time in fm/c
	         G4double time = currentTime + decayTime;
	         if(time <= maximumTime) {
	           registerAvatar(new DecayAvatar((*i), time, theNucleus));
	         }
	       }
        if((*i)->isOmega()) {
          G4double decayTimeOmega = PionResonanceDecayChannel::computeDecayTime((*i));
          G4double timeOmega = currentTime + decayTimeOmega;
          if(timeOmega <= maximumTime) {
            registerAvatar(new DecayAvatar((*i), timeOmega, theNucleus));
          }
        }
      }
    }

    G4INCL::IAvatar* StandardPropagationModel::propagate(FinalState const * const fs)
    {
      if(fs) {
        // We update only the information related to particles that were updated
        // by the previous avatar.
#ifdef INCL_REGENERATE_AVATARS
#warning "The INCL_REGENERATE_AVATARS code has not been tested in a while. Use it at your peril."
        if(!fs->getModifiedParticles().empty() || !fs->getEnteringParticles().empty() || !fs->getCreatedParticles().empty()) {
          // Regenerates the entire avatar list, skipping collisions between
          // updated particles
          theNucleus->getStore()->clearAvatars();
          theNucleus->getStore()->initialiseParticleAvatarConnections();
          generateAllAvatarsExceptUpdated(fs);
        }
#else
        ParticleList const &updatedParticles = fs->getModifiedParticles();
        if(fs->getValidity()==PauliBlockedFS) {
          // This final state might represents the outcome of a Pauli-blocked delta
          // decay
// assert(updatedParticles.empty() || (updatedParticles.size()==1 && updatedParticles.front()->isResonance()));
// assert(fs->getEnteringParticles().empty() && fs->getCreatedParticles().empty() && fs->getOutgoingParticles().empty() && fs->getDestroyedParticles().empty());
          generateDecays(updatedParticles);
        } else {
          ParticleList const &entering = fs->getEnteringParticles();
          generateDecays(updatedParticles);
          generateDecays(entering);

          ParticleList const &created = fs->getCreatedParticles();
          if(created.empty() && entering.empty())
            updateAvatars(updatedParticles);
          else {
            ParticleList updatedParticlesCopy = updatedParticles;
            updatedParticlesCopy.insert(updatedParticlesCopy.end(), entering.begin(), entering.end());
            updatedParticlesCopy.insert(updatedParticlesCopy.end(), created.begin(), created.end());
            updateAvatars(updatedParticlesCopy);
          }
        }
#endif
      }

      G4INCL::IAvatar *theAvatar = theNucleus->getStore()->findSmallestTime();
      if(theAvatar == 0) return 0; // Avatar list is empty
      //      theAvatar->dispose();

      if(theAvatar->getTime() < currentTime) {
        INCL_ERROR("Avatar time = " << theAvatar->getTime() << ", currentTime = " << currentTime << '\n');
        return 0;
      } else if(theAvatar->getTime() > currentTime) {
        theNucleus->getStore()->timeStep(theAvatar->getTime() - currentTime);

        currentTime = theAvatar->getTime();
        theNucleus->getStore()->getBook().setCurrentTime(currentTime);
      }

      return theAvatar;
    }
}
