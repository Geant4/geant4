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
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/*
 * StandardPropagationModel.cpp
 *
 *  Created on: 4 juin 2009
 *      Author: Pekka Kaitaniemi
 */

#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLSurfaceAvatar.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLDecayAvatar.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLRandom.hh"
#include <iostream>
#include "G4INCLLogger.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLCoulombDistortion.hh"
#include "G4INCLDeltaDecayChannel.hh"

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

    G4bool StandardPropagationModel::shootProjectile(G4INCL::Particle *p, G4double impactParameter) {
      firstAvatar = true;
      currentTime = 0.0;

      G4double temfin = 0.0;
      if( p->isNucleon() )
        temfin = 29.8 * std::pow(theNucleus->getA(), 0.16);
      else {
        const G4double tlab = p->getEnergy() - p->getMass();
        temfin = 30.18 * std::pow(theNucleus->getA(), 0.17*(1.0 - 5.7E-5*tlab));
      }

      maximumTime = temfin;

      // If Coulomb is activated, do not process events with impact
      // parameter larger than the maximum impact parameter, taking G4into
      // account Coulomb distortion.
      if(impactParameter>CoulombDistortion::maxImpactParameter(p,theNucleus))
        return false;

      const G4double tbid = Random::shoot() * Math::twoPi;
      ThreeVector position(impactParameter * std::cos(tbid),
          impactParameter * std::sin(tbid),
          -1.E3);
      p->setPosition(position);

      theNucleus->setIncomingAngularMomentum(p->getAngularMomentum());
      theNucleus->setIncomingMomentum(p->getMomentum());
      theNucleus->setInitialEnergy(p->getEnergy() + ParticleTable::getMass(theNucleus->getA(),theNucleus->getZ()));

      CoulombDistortion::bringToSurface(p, theNucleus);

      theNucleus->getStore()->addIncomingParticle(p); // puts the particle in the waiting list
      theNucleus->particleEnters(p); // removes the particle from the waiting list and changes its kinetic energy
      theNucleus->insertParticipant(p);
      return true;
    }

    G4bool StandardPropagationModel::shootProjectile(G4INCL::Nucleus * /* p */, G4double /* impactParameter */) {
      firstAvatar = true;
      currentTime = 0.0;
      return true;
    }

    G4double StandardPropagationModel::getStoppingTime() {
      return maximumTime;
    }

    void StandardPropagationModel::setStoppingTime(G4double time) {
      if(time > 0.0) {
        maximumTime = time;
      } else {
        ERROR("new stopping time is smaller than 0!" << std::endl);
      }
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
      if(!p1->isParticipant() && !p2->isParticipant()) return NULL;

      // Is it a pi-resonance collision (we don't treat them)?
      if((p1->isResonance() && p2->isPion()) || (p1->isPion() && p2->isResonance()))
        return NULL;

      // 2N < Tf
      if(
          (p1->isNucleon() && p1->getKineticEnergy()<theNucleus->getPotential()->getFermiEnergy(p1->getType())) &&
          (p2->isNucleon() && p2->getKineticEnergy()<theNucleus->getPotential()->getFermiEnergy(p2->getType()))
          )
        return NULL;

      // Is the CM energy > cutNN? (no cutNN on the first collision)
      if(theNucleus->getStore()->getBook()->getAcceptedCollisions()>0
          && p1->isNucleon() && p2->isNucleon()
          && KinematicsUtils::squareTotalEnergyInCM(p1,p2) < BinaryCollisionAvatar::cutNNSquared) return NULL;

      // Will the avatar take place between now and the end of the cascade?
      G4double minDistOfApproachSquared = 0.0;
      G4double t = getTime(p1, p2, &minDistOfApproachSquared);
      if(t>maximumTime || t<currentTime) return NULL;

      // Local energy. Jump through some hoops to calculate the cross section
      // at the collision poG4int, and clean up after yourself afterwards.
      ThreeVector mom1, mom2, pos1, pos2;
      G4double energy1 = 0.0, energy2 = 0.0;
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

      if(p1HasLocalEnergy) {
        mom1 = p1->getMomentum();
        pos1 = p1->getPosition();
        energy1 = p1->getEnergy();
        p1->propagate(t - currentTime);
        if(p1->getPosition().mag() > theNucleus->getSurfaceRadius(p1)) {
          p1->setPosition(pos1);
          p1->setMomentum(mom1);
          p1->setEnergy(energy1);
          return NULL;
        }
        KinematicsUtils::transformToLocalEnergyFrame(theNucleus, p1);
      }
      if(p2HasLocalEnergy) {
        energy2 = p2->getEnergy();
        mom2 = p2->getMomentum();
        pos2 = p2->getPosition();
        p2->propagate(t - currentTime);
        if(p2->getPosition().mag() > theNucleus->getSurfaceRadius(p2)) {
          p2->setPosition(pos2);
          p2->setMomentum(mom2);
          p2->setEnergy(energy2);
          if(p1HasLocalEnergy) {
            p1->setPosition(pos1);
            p1->setMomentum(mom1);
            p1->setEnergy(energy1);
          }
          return NULL;
        }
        KinematicsUtils::transformToLocalEnergyFrame(theNucleus, p2);
      }

      // Compute the total cross section
      const G4double totalCrossSection = CrossSections::total(p1, p2);

      // Restore particles to their state before the local-energy tweak
      if(p1HasLocalEnergy) {
        p1->setPosition(pos1);
        p1->setMomentum(mom1);
        p1->setEnergy(energy1);
      }
      if(p2HasLocalEnergy) {
        p2->setPosition(pos2);
        p2->setMomentum(mom2);
        p2->setEnergy(energy2);
      }

      // Do the particles come close enough to each other?
      if(Math::tenPi*minDistOfApproachSquared > totalCrossSection) return NULL;

      // Warn if the two collision partners are the same particle
      if(p1->getID() == p2->getID()) {
        ERROR("At BinaryCollisonAvatar generation, ID1 (" << p1->getID()
          << ") == ID2 (" << p2->getID() <<")." << std::endl);
      }

      // Return a new avatar, then!
      return new G4INCL::BinaryCollisionAvatar(t, totalCrossSection, theNucleus, p1, p2);
    }

    G4double StandardPropagationModel::getReflectionTime(G4INCL::Particle const * const aParticle) {
      G4double time = 0.0;
      const G4double T2 = aParticle->getMomentum().mag2();
      const G4double T4 = aParticle->getPosition().mag2();
      
      const G4double r = theNucleus->getSurfaceRadius(aParticle);
      const G4double T1 = aParticle->getMomentum().dot(aParticle->getPosition());
      const G4double T3 = T1/T2;
      const G4double T5 = T3*T3 + (r*r-T4)/T2;
      if(T5 < 0.0) {
        ERROR("Imaginary reflection time! Delta = " << T5 << " for particle: " << std::endl
          << aParticle->prG4int());
        time = 10000.0;
      } else {
        time = currentTime + (-T3 + std::sqrt(T5)) * aParticle->getEnergy();
      }
      return time;
    }

    G4double StandardPropagationModel::getTime(G4INCL::Particle const * const particleA,
					     G4INCL::Particle const * const particleB, G4double *minDistOfApproach) const
    {
      G4double time;
      G4INCL::ThreeVector t13 = particleA->getMomentum()/particleA->getEnergy();
      t13 -= particleB->getMomentum()/particleB->getEnergy();
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

    void StandardPropagationModel::checkCollisions(const ParticleList &participants,
        const ParticleList &particles)
    {
      G4int iind = 1, jind = 1;
      for(ParticleIter i = participants.begin(); i != participants.end(); ++i) {
        jind = 0;
        for(ParticleIter j = particles.begin(); j != particles.end(); ++j) {
          if( (*j)->isParticipant() )
          {
            ++jind;
            if(jind >= iind) continue;
          }
          if((*i)->getID() == (*j)->getID()) continue; // Do not process the collision of a particle with itself

          registerAvatar(generateBinaryCollisionAvatar(*i,*j));
        }
        ++iind;
      }
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
      ParticleList p = theNucleus->getStore()->getParticles();
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
      if(firstAvatar) { // When we propagate particles for the first time we create the full list of avatars.
	generateAllAvatars();
/*        if(!theNucleus->getStore()->containsCollisions()) {
          theNucleus->forceTransparent();
          return NULL;
        }*/
        firstAvatar = false;
      } else { // For subsequent avatars we update only the
        // information related to particles that were updated
        // by the previous avatar.
#ifdef INCL_REGENERATE_AVATARS
#warning "The INCL_REGENERATE_AVATARS code has not been tested in a while. Use it at your peril."
        // Regenerates the entire avatar list, skipping collisions between
        // updated particles
        if(theNucleus->getUpdatedParticles().size()!=0 || theNucleus->getCreatedParticles().size()!=0) {
          theNucleus->getStore()->clearAvatars();
          theNucleus->getStore()->initialiseParticleAvatarConnections();
          generateAllAvatars(true);
        }
#else
        // Deltas are created by transforming nucleon G4into a delta for
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
      }

      G4INCL::IAvatar *theAvatar = theNucleus->getStore()->findSmallestTime();
      if(theAvatar == 0) return 0; // Avatar list is empty
      //      theAvatar->dispose();

      theNucleus->getStore()->timeStep(theAvatar->getTime() - currentTime);

      if(theAvatar->getTime() <= currentTime) {
        ERROR("Avatar time = " << theAvatar->getTime() << ", currentTime = " << currentTime << std::endl);
        return 0;
      } else {
        currentTime = theAvatar->getTime();
        theNucleus->getStore()->getBook()->setCurrentTime(currentTime);
      }

      return theAvatar;
    }
}
