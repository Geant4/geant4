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
 * G4INCLBinaryCollisionAvatar.cc
 *
 *  \date Jun 5, 2009
 * \author Pekka Kaitaniemi
 */

#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLElasticChannel.hh"
#include "G4INCLRecombinationChannel.hh"
#include "G4INCLDeltaProductionChannel.hh"
#include "G4INCLNNToMultiPionsChannel.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLRandom.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLPauliBlocking.hh"
#include "G4INCLPiNElasticChannel.hh"
#include "G4INCLPiNToDeltaChannel.hh"
#include "G4INCLPiNToMultiPionsChannel.hh"
#include "G4INCLStore.hh"
#include "G4INCLBook.hh"
#include "G4INCLLogger.hh"
#include <string>
#include <sstream>
// #include <cassert>

namespace G4INCL {

  // WARNING: if you update the default cutNN value, make sure you update the
  // cutNNSquared variable, too.
  G4ThreadLocal G4double BinaryCollisionAvatar::cutNN = 1910.0;
  G4ThreadLocal G4double BinaryCollisionAvatar::cutNNSquared = 3648100.0; // 1910.0 * 1910.0

  BinaryCollisionAvatar::BinaryCollisionAvatar(G4double time, G4double crossSection,
      G4INCL::Nucleus *n, G4INCL::Particle *p1, G4INCL::Particle *p2)
    : InteractionAvatar(time, n, p1, p2), theCrossSection(crossSection),
    isParticle1Spectator(false),
    isParticle2Spectator(false),
    isElastic(false)
  {
    setType(CollisionAvatarType);
  }

  BinaryCollisionAvatar::~BinaryCollisionAvatar() {
  }

  G4INCL::IChannel* BinaryCollisionAvatar::getChannel() {
    // We already check cutNN at avatar creation time, but we have to check it
    // again here. For composite projectiles, we might have created independent
    // avatars with no cutNN before any collision took place.
    if(particle1->isNucleon()
        && particle2->isNucleon()
        && theNucleus->getStore()->getBook().getAcceptedCollisions()!=0) {
      const G4double energyCM2 = KinematicsUtils::squareTotalEnergyInCM(particle1, particle2);
      // Below a certain cut value we don't do anything:
      if(energyCM2 < cutNNSquared) {
        INCL_DEBUG("CM energy = sqrt(" << energyCM2 << ") MeV < sqrt(" << cutNNSquared
            << ") MeV = cutNN" << "; returning a NULL channel" << '\n');
        InteractionAvatar::restoreParticles();
        return NULL;
      }
    }

    /** Check again the distance of approach. In order for the avatar to be
     * realised, we have to perform a check in the CM system. We define a
     * distance four-vector as
     * \f[ (0, \Delta\vec{x}), \f]
     * where \f$\Delta\vec{x}\f$ is the distance vector of the particles at
     * their minimum distance of approach (i.e. at the avatar time). By
     * boosting this four-vector to the CM frame of the two particles and we
     * obtain a new four vector
     * \f[ (\Delta t', \Delta\vec{x}'), \f]
     * with a non-zero time component (the collision happens simultaneously for
     * the two particles in the lab system, but not in the CM system). In order
     * for the avatar to be realised, we require that
     * \f[ |\Delta\vec{x}'| \leq \sqrt{\sigma/\pi}.\f]
     * Note that \f$|\Delta\vec{x}'|\leq|\Delta\vec{x}|\f$; thus, the condition
     * above is more restrictive than the check that we perform in
     * G4INCL::Propagation::StandardPropagationModel::generateBinaryCollisionAvatar.
     * In other words, the avatar generation cannot miss any physical collision
     * avatars.
     */
    ThreeVector minimumDistance = particle1->getPosition();
    minimumDistance -= particle2->getPosition();
    const G4double betaDotX = boostVector.dot(minimumDistance);
    const G4double minDist = Math::tenPi*(minimumDistance.mag2() + betaDotX*betaDotX / (1.-boostVector.mag2()));
    if(minDist > theCrossSection) {
      INCL_DEBUG("CM distance of approach is too small: " << minDist << ">" <<
        theCrossSection <<"; returning a NULL channel" << '\n');
      InteractionAvatar::restoreParticles();
      return NULL;
    }

//// NN
    if(particle1->isNucleon() && particle2->isNucleon()) {
      const G4double elasticCX = CrossSections::elastic(particle1, particle2);
      const G4double deltaProductionCX = CrossSections::NNToNDelta(particle1, particle2);
      const G4double onePiProductionCX = CrossSections::NNToxPiNN(1,particle1, particle2);
      const G4double twoPiProductionCX = CrossSections::NNToxPiNN(2,particle1, particle2);
      const G4double threePiProductionCX = CrossSections::NNToxPiNN(3,particle1, particle2);
      const G4double fourPiProductionCX = CrossSections::NNToxPiNN(4,particle1, particle2);
      const G4double totCX=CrossSections::total(particle1, particle2);

      const G4double rChannel=Random::shoot() * totCX;

      if(elasticCX > rChannel) {
        // Elastic NN channel
        isElastic = true;
        INCL_DEBUG("NN interaction: elastic channel chosen" << '\n');
        return new ElasticChannel(particle1, particle2);
      } else if((elasticCX + deltaProductionCX) > rChannel) {
        isElastic = false;
        // NN -> N Delta channel is chosen
        INCL_DEBUG("NN interaction: Delta channel chosen" << '\n');
        return new DeltaProductionChannel(particle1, particle2);
      } else if(elasticCX + deltaProductionCX + onePiProductionCX > rChannel) {
        isElastic = false;
        // NN -> PiNN channel is chosen
        INCL_DEBUG("NN interaction: one Pion channel chosen" << '\n');
        return new NNToMultiPionsChannel(1,particle1, particle2);
      } else if(elasticCX + deltaProductionCX + onePiProductionCX + twoPiProductionCX > rChannel) {
        isElastic = false;
        // NN -> 2PiNN channel is chosen
        INCL_DEBUG("NN interaction: two Pions channel chosen" << '\n');
        return new NNToMultiPionsChannel(2,particle1, particle2);
      } else if(elasticCX + deltaProductionCX + onePiProductionCX + twoPiProductionCX + threePiProductionCX > rChannel) {
        isElastic = false;
        // NN -> 3PiNN channel is chosen
        INCL_DEBUG("NN interaction: three Pions channel chosen" << '\n');
        return new NNToMultiPionsChannel(3,particle1, particle2);
      } else if (elasticCX + deltaProductionCX + onePiProductionCX + twoPiProductionCX + threePiProductionCX + fourPiProductionCX > rChannel) {
        isElastic = false;
        // NN -> 4PiNN channel is chosen
        INCL_DEBUG("NN interaction: four Pions channel chosen" << '\n');
        return new NNToMultiPionsChannel(4,particle1, particle2);
      } else {
        INCL_WARN("inconsistency with NN Cross Sections: returning an elastic channel" << '\n');
        isElastic = true;
        return new ElasticChannel(particle1, particle2);
      }

//// NDelta
    } else if((particle1->isNucleon() && particle2->isDelta()) ||
	      (particle1->isDelta() && particle2->isNucleon())) {
          G4double elasticCX = CrossSections::elastic(particle1, particle2);
          G4double recombinationCX = CrossSections::NDeltaToNN(particle1, particle2);

          if(elasticCX/(elasticCX + recombinationCX) < Random::shoot()) {
            isElastic = false;
          } else
            isElastic = true;

          if(isElastic) {
// Elastic N Delta channel
              INCL_DEBUG("NDelta interaction: elastic channel chosen" << '\n');
              return new ElasticChannel(particle1, particle2);
          } else { // Recombination
// N Delta -> NN channel is chosen
              INCL_DEBUG("NDelta interaction: recombination channel chosen" << '\n');
              return new RecombinationChannel(particle1, particle2);
            }

//// DeltaDelta
    } else if(particle1->isDelta() && particle2->isDelta()) {
        isElastic = true;
        INCL_DEBUG("DeltaDelta interaction: elastic channel chosen" << '\n');
        return new ElasticChannel(particle1, particle2);

//// PiN
    } else if(isPiN) {
      const G4double elasticCX = CrossSections::elastic(particle1, particle2);
      const G4double deltaProductionCX = CrossSections::piNToDelta(particle1, particle2);
      const G4double onePiProductionCX = CrossSections::piNToxPiN(2,particle1, particle2);
      const G4double twoPiProductionCX = CrossSections::piNToxPiN(3,particle1, particle2);
      const G4double threePiProductionCX = CrossSections::piNToxPiN(4,particle1, particle2);
      const G4double totCX=CrossSections::total(particle1, particle2);
// assert(std::fabs(totCX-elasticCX-deltaProductionCX-onePiProductionCX-twoPiProductionCX-threePiProductionCX)<1.);

      const G4double rChannel=Random::shoot() * totCX;

      if(elasticCX > rChannel) {
        // Elastic PiN channel
        isElastic = true;
        INCL_DEBUG("PiN interaction: elastic channel chosen" << '\n');
        return new PiNElasticChannel(particle1, particle2);
      } else if(elasticCX + deltaProductionCX > rChannel) {
        isElastic = false;
        // PiN -> Delta channel is chosen
        INCL_DEBUG("PiN interaction: Delta channel chosen" << '\n');
        return new PiNToDeltaChannel(particle1, particle2);
      } else if(elasticCX + deltaProductionCX + onePiProductionCX > rChannel) {
        isElastic = false;
        // PiN -> PiNPi channel is chosen
        INCL_DEBUG("PiN interaction: one Pion channel chosen" << '\n');
        return new PiNToMultiPionsChannel(2,particle1, particle2);
      } else if(elasticCX + deltaProductionCX + onePiProductionCX + twoPiProductionCX > rChannel) {
        isElastic = false;
        // PiN -> PiN2Pi channel is chosen
        INCL_DEBUG("PiN interaction: two Pions channel chosen" << '\n');
        return new PiNToMultiPionsChannel(3,particle1, particle2);
      } else if(elasticCX + deltaProductionCX + onePiProductionCX + twoPiProductionCX + threePiProductionCX > rChannel) {
        isElastic = false;
        // PiN -> PiN3Pi channel is chosen
        INCL_DEBUG("PiN interaction: three Pions channel chosen" << '\n');
        return new PiNToMultiPionsChannel(4,particle1, particle2);
      } else {
        INCL_WARN("inconsistency with PiN Cross Sections: returning an elastic channel" << '\n');
        isElastic = true;
        return new PiNElasticChannel(particle1, particle2);
      }
    } else {
      INCL_DEBUG("BinaryCollisionAvatar can only handle nucleons (for the moment)."
	      << '\n'
	      << particle1->print()
	      << '\n'
	      << particle2->print()
	      << '\n');
      InteractionAvatar::restoreParticles();
      return NULL;
    }
  }

  void BinaryCollisionAvatar::preInteraction() {
    isParticle1Spectator = particle1->isTargetSpectator();
    isParticle2Spectator = particle2->isTargetSpectator();
    InteractionAvatar::preInteraction();
  }

  void BinaryCollisionAvatar::postInteraction(FinalState *fs) {
    // Call the postInteraction method of the parent class
    // (provides Pauli blocking and enforces energy conservation)
    InteractionAvatar::postInteraction(fs);

    switch(fs->getValidity()) {
      case PauliBlockedFS:
        theNucleus->getStore()->getBook().incrementBlockedCollisions();
        break;
      case NoEnergyConservationFS:
      case ParticleBelowFermiFS:
      case ParticleBelowZeroFS:
        break;
      case ValidFS:
        Book &theBook = theNucleus->getStore()->getBook();
        theBook.incrementAcceptedCollisions();
        if(theBook.getAcceptedCollisions() == 1) {
          // Store time and cross section of the first collision
          G4double t = theBook.getCurrentTime();
          theBook.setFirstCollisionTime(t);
          theBook.setFirstCollisionXSec(oldXSec);

          // Store position and momentum of the spectator on the first
          // collision
          if((isParticle1Spectator && isParticle2Spectator) || (!isParticle1Spectator && !isParticle2Spectator)) {
            INCL_ERROR("First collision must be within a target spectator and a non-target spectator");
          }
          if(isParticle1Spectator) {
            theBook.setFirstCollisionSpectatorPosition(backupParticle1->getPosition().mag());
            theBook.setFirstCollisionSpectatorMomentum(backupParticle1->getMomentum().mag());
          } else {
            theBook.setFirstCollisionSpectatorPosition(backupParticle2->getPosition().mag());
            theBook.setFirstCollisionSpectatorMomentum(backupParticle2->getMomentum().mag());
          }

          // Store the elasticity of the first collision
          theBook.setFirstCollisionIsElastic(isElastic);
        }
    }
    return;
  }

  std::string BinaryCollisionAvatar::dump() const {
    std::stringstream ss;
    ss << "(avatar " << theTime <<" 'nn-collision" << '\n'
      << "(list " << '\n'
      << particle1->dump()
      << particle2->dump()
      << "))" << '\n';
    return ss.str();
  }

}
