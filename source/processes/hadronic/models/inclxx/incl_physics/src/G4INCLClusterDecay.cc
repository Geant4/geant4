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

/** \file G4INCLClusterDecay.cc
 * \brief Static class for carrying out cluster decays
 *
 * Created on: 6th July 2011
 *     Author: Davide Mancusi
 */

#include "G4INCLClusterDecay.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLRandom.hh"
//#include <cassert>

namespace G4INCL {

  ParticleList ClusterDecay::decay(Cluster * const c) {
    ParticleList decayProducts;
    recursiveDecay(c, &decayProducts);
    return decayProducts;
  }

  void ClusterDecay::recursiveDecay(Cluster * const c, ParticleList *decayProducts) {
    const G4int Z = c->getZ();
    const G4int A = c->getA();

    ParticleTable::ClusterDecayType theDecayMode = ParticleTable::clusterDecayMode[Z][A];

    switch(theDecayMode) {
      default:
        ERROR("Unrecognized cluster-decay mode: " << theDecayMode << std::endl
            << c->prG4int());
      case ParticleTable::StableCluster:
        // For stable clusters, just return
        return;
        break;
      case ParticleTable::ProtonDecay:
      case ParticleTable::NeutronDecay:
      case ParticleTable::AlphaDecay:
        // Two-body decays
        twoBodyDecay(c, theDecayMode, decayProducts);
        break;
      case ParticleTable::TwoProtonDecay:
      case ParticleTable::TwoNeutronDecay:
        // Three-body decays
        threeBodyDecay(c, theDecayMode, decayProducts);
        break;
    }

    // Calls itself recursively in case the produced remnant is still unstable.
    // Sneaky, isn't it.
    recursiveDecay(c,decayProducts);
  }

  void ClusterDecay::twoBodyDecay(Cluster * const c, ParticleTable::ClusterDecayType theDecayMode, ParticleList *decayProducts) {
    Particle *decayParticle = 0;
    const ThreeVector mom(0.0, 0.0, 0.0);
    const ThreeVector pos = c->getPosition();

    // Create the emitted particle
    switch(theDecayMode) {
      case ParticleTable::ProtonDecay:
        decayParticle = new Particle(Proton, mom, pos);
        break;
      case ParticleTable::NeutronDecay:
        decayParticle = new Particle(Neutron, mom, pos);
        break;
      case ParticleTable::AlphaDecay:
        decayParticle = new Cluster(2,4);
        break;
      default:
        ERROR("Unrecognized cluster-decay mode in two-body decay: " << theDecayMode << std::endl
            << c->prG4int());
        return;
    }
    decayParticle->makeParticipant();
    decayParticle->setNumberOfDecays(1);
    decayParticle->setPosition(c->getPosition());
    decayParticle->setEmissionTime(c->getEmissionTime());

    // Save some variables of the mother cluster
    const G4double motherMass = c->getMass();
    const ThreeVector velocity = -c->boostVector();

    // Characteristics of the daughter particle
    const G4int daughterZ = c->getZ() - decayParticle->getZ();
    const G4int daughterA = c->getA() - decayParticle->getA();
    const G4double daughterMass = ParticleTable::getMass(daughterA,daughterZ);

    // The mother cluster becomes the daughter
    c->setZ(daughterZ);
    c->setA(daughterA);
    c->setMass(daughterMass);

    const G4double decayMass = decayParticle->getMass();
    // assert(motherMass > daughterMass + decayMass); // Q-value should be >0

    // Decay kinematics in the mother rest frame
    const G4double pCM = KinematicsUtils::momentumInCM(motherMass, daughterMass, decayMass);
    const ThreeVector momentum = Random::normVector(pCM);
    c->setMomentum(momentum);
    c->adjustEnergyFromMomentum();
    decayParticle->setMomentum(-momentum);
    decayParticle->adjustEnergyFromMomentum();

    // Boost to the lab frame
    decayParticle->boost(velocity);
    c->boost(velocity);

    // Add the decay particle to the list of decay products
    decayProducts->push_back(decayParticle);
  }

  void ClusterDecay::threeBodyDecay(Cluster * const c, ParticleTable::ClusterDecayType theDecayMode, ParticleList *decayProducts) {
    Particle *decayParticle1 = 0;
    Particle *decayParticle2 = 0;
    const ThreeVector mom(0.0, 0.0, 0.0);
    const ThreeVector pos = c->getPosition();

    // Create the emitted particles
    switch(theDecayMode) {
      case ParticleTable::TwoProtonDecay:
        decayParticle1 = new Particle(Proton, mom, pos);
        decayParticle2 = new Particle(Proton, mom, pos);
        break;
      case ParticleTable::TwoNeutronDecay:
        decayParticle1 = new Particle(Neutron, mom, pos);
        decayParticle2 = new Particle(Neutron, mom, pos);
        break;
      default:
        ERROR("Unrecognized cluster-decay mode in three-body decay: " << theDecayMode << std::endl
            << c->prG4int());
        return;
    }
    decayParticle1->makeParticipant();
    decayParticle2->makeParticipant();
    decayParticle1->setNumberOfDecays(1);
    decayParticle2->setNumberOfDecays(1);

    // Save some variables of the mother cluster
    const G4double motherMass = c->getMass();
    const ThreeVector velocity = -c->boostVector();

    // Masses and charges of the daughter particle and of the decay products
    const G4int decayZ1 = decayParticle1->getZ();
    const G4int decayA1 = decayParticle1->getA();
    const G4int decayZ2 = decayParticle2->getZ();
    const G4int decayA2 = decayParticle2->getA();
    const G4int decayZ = decayZ1 + decayZ2;
    const G4int decayA = decayA1 + decayA2;
    const G4int daughterZ = c->getZ() - decayZ;
    const G4int daughterA = c->getA() - decayA;
    const G4double decayMass1 = decayParticle1->getMass();
    const G4double decayMass2 = decayParticle2->getMass();
    const G4double daughterMass = ParticleTable::getMass(daughterA,daughterZ);

    // Q-values
    const G4double qValue = motherMass - daughterMass - decayMass1 - decayMass2;
    // assert(qValue > 0.); // Q-value should be >0
    const G4double qValueB = qValue * Random::shoot();

    // The decay particles behave as if they had more mass until the second
    // decay
    const G4double decayMass = decayMass1 + decayMass2 + qValueB;

    /* Stage A: mother --> daughter + (decay1+decay2) */

    // The mother cluster becomes the daughter
    c->setZ(daughterZ);
    c->setA(daughterA);
    c->setMass(daughterMass);

    // Decay kinematics in the mother rest frame
    const G4double pCMA = KinematicsUtils::momentumInCM(motherMass, daughterMass, decayMass);
    const ThreeVector momentumA = Random::normVector(pCMA);
    c->setMomentum(momentumA);
    c->adjustEnergyFromMomentum();
    const ThreeVector decayBoostVector = momentumA/std::sqrt(decayMass*decayMass + momentumA.mag2());

    /* Stage B: (decay1+decay2) --> decay1 + decay2 */

    // Decay kinematics in the (decay1+decay2) rest frame
    const G4double pCMB = KinematicsUtils::momentumInCM(decayMass, decayMass1, decayMass2);
    const ThreeVector momentumB = Random::normVector(pCMB);
    decayParticle1->setMomentum(momentumB);
    decayParticle2->setMomentum(-momentumB);
    decayParticle1->adjustEnergyFromMomentum();
    decayParticle2->adjustEnergyFromMomentum();

    // Boost decay1 and decay2 to the Stage-A decay frame
    decayParticle1->boost(decayBoostVector);
    decayParticle2->boost(decayBoostVector);

    // Boost all particles to the lab frame
    decayParticle1->boost(velocity);
    decayParticle2->boost(velocity);
    c->boost(velocity);

    // Add the decay particles to the list of decay products
    decayProducts->push_back(decayParticle1);
    decayProducts->push_back(decayParticle2);
  }

}

