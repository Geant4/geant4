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
 * G4INCLBinaryCollisionAvatar.cc
 *
 *  \date Jun 5, 2009
 * \author Pekka Kaitaniemi
 */

#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLElasticChannel.hh"
#include "G4INCLRecombinationChannel.hh"
#include "G4INCLDeltaProductionChannel.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLRandom.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLPauliBlocking.hh"
#include "G4INCLPionNucleonChannel.hh"
#include "G4INCLStore.hh"
#include "G4INCLBook.hh"
#include "G4INCLLogger.hh"
#include <string>
#include <sstream>
// #include <cassert>

namespace G4INCL {

  const G4double BinaryCollisionAvatar::cutNN = 1910;
  const G4double BinaryCollisionAvatar::cutNNSquared = cutNN*cutNN;

  BinaryCollisionAvatar::BinaryCollisionAvatar(G4double time, G4double crossSection,
      G4INCL::Nucleus *n, G4INCL::Particle *p1, G4INCL::Particle *p2)
    : InteractionAvatar(time, n, p1, p2), theCrossSection(crossSection)
  {
    setType(CollisionAvatarType);
  }

  BinaryCollisionAvatar::~BinaryCollisionAvatar() {
  }

  G4INCL::IChannel* BinaryCollisionAvatar::getChannel() const {
    // We already check cutNN at avatar creation time, but we have to check it
    // again here. For composite projectiles, we might have created independent
    // avatars with no cutNN before any collision took place.
    if(particle1->isNucleon()
        && particle2->isNucleon()
        && theNucleus->getStore()->getBook()->getAcceptedCollisions()!=0) {
      const G4double energyCM2 = KinematicsUtils::squareTotalEnergyInCM(particle1, particle2);
      // Below a certain cut value we don't do anything:
      if(energyCM2 < cutNNSquared) {
        DEBUG("CM energy = sqrt(" << energyCM2 << ") MeV < sqrt(" << cutNNSquared
            << ") MeV = cutNN" << "; returning a NULL channel" << std::endl);
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
      DEBUG("CM distance of approach is too small: " << minDist << ">" <<
        theCrossSection <<"; returning a NULL channel" << std::endl);
      InteractionAvatar::restoreParticles();
      return NULL;
    }

    if(particle1->isNucleon() && particle2->isNucleon()) { // NN->NN
      G4double elasticCX = CrossSections::elastic(particle1,
          particle2);
      G4double deltaProductionCX = CrossSections::deltaProduction(particle1,
          particle2);

      G4bool isElastic = true;
      if(elasticCX/(elasticCX + deltaProductionCX) < Random::shoot()) {
        // NN -> N Delta channel is chosen
        isElastic = false;
      }

      if(isElastic) { // Elastic NN channel
        DEBUG("NN interaction: elastic channel chosen" << std::endl);
        return new ElasticChannel(theNucleus, particle1, particle2);
      } else { // Delta production
        // Inelastic NN channel
        DEBUG("NN interaction: inelastic channel chosen" << std::endl);
        return new DeltaProductionChannel(particle1, particle2, theNucleus);
      }
    } else if((particle1->isNucleon() && particle2->isDelta()) ||
	      (particle1->isDelta() && particle2->isNucleon())) {
      G4double elasticCX = CrossSections::elastic(particle1,
          particle2);
      G4double recombinationCX = CrossSections::recombination(particle1,
          particle2);

      G4bool isElastic = true;
      if(elasticCX/(elasticCX + recombinationCX) < Random::shoot()) {
        // N Delta -> NN channel is chosen
        isElastic = false;
      }

      if(isElastic) { // Elastic N Delta channel
        DEBUG("NDelta interaction: elastic channel chosen" << std::endl);
        return new ElasticChannel(theNucleus, particle1, particle2);
      } else { // Recombination
        DEBUG("NDelta interaction: recombination channel chosen" << std::endl);
        return new RecombinationChannel(theNucleus, particle1, particle2);
      }
    } else if(particle1->isDelta() && particle2->isDelta()) {
        DEBUG("DeltaDelta interaction: elastic channel chosen" << std::endl);
        return new ElasticChannel(theNucleus, particle1, particle2);
    } else if((particle1->isNucleon() && particle2->isPion()) ||
	      (particle1->isPion() && particle2->isNucleon())) {
      return new PionNucleonChannel(particle1, particle2, theNucleus, shouldUseLocalEnergy());
    } else {
      DEBUG("BinaryCollisionAvatar can only handle nucleons (for the moment)."
	      << std::endl
	      << particle1->print()
	      << std::endl
	      << particle2->print()
	      << std::endl);
      InteractionAvatar::restoreParticles();
      return NULL;
    }
  }

  void BinaryCollisionAvatar::preInteraction() {
    InteractionAvatar::preInteraction();
  }

  FinalState *BinaryCollisionAvatar::postInteraction(FinalState *fs) {
    // Call the postInteraction method of the parent class
    // (provides Pauli blocking and enforces energy conservation)
    fs = InteractionAvatar::postInteraction(fs);

    switch(fs->getValidity()) {
      case PauliBlockedFS:
        theNucleus->getStore()->getBook()->incrementBlockedCollisions();
        break;
      case NoEnergyConservationFS:
      case ParticleBelowFermiFS:
      case ParticleBelowZeroFS:
        break;
      case ValidFS:
        theNucleus->getStore()->getBook()->incrementAcceptedCollisions();
        if(theNucleus->getStore()->getBook()->getAcceptedCollisions() == 1) {
          G4double t = theNucleus->getStore()->getBook()->getCurrentTime();
          theNucleus->getStore()->getBook()->setFirstCollisionTime(t);
          theNucleus->getStore()->getBook()->setFirstCollisionXSec(oldXSec);
        }
    }
    return fs;
  }

  std::string BinaryCollisionAvatar::dump() const {
    std::stringstream ss;
    ss << "(avatar " << theTime <<" 'nn-collision" << std::endl
      << "(list " << std::endl
      << particle1->dump()
      << particle2->dump()
      << "))" << std::endl;
    return ss.str();
  }

}
