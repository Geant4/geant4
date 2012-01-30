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
// INCL++ revision: v5.1_rc1
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/*
 * G4INCLNucleus.cc
 *
 *  Created on: Jun 5, 2009
 *      Author: Pekka Kaitaniemi
 */

#ifndef G4INCLNucleus_hh
#define G4INCLNucleus_hh 1

#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include "G4INCLParticle.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLNuclearDensityFactory.hh"
#include "G4INCLNuclearPotentialConstant.hh"
#include "G4INCLNuclearPotentialIsospin.hh"
#include "G4INCLNuclearPotentialEnergyIsospin.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLDecayAvatar.hh"
#include "G4INCLCluster.hh"
#include "G4INCLClusterDecay.hh"
#include <iterator>
#include <cstdlib>
#include <sstream>
// #include <cassert>

namespace G4INCL {

  Nucleus::Nucleus(G4int mass, G4int charge, Config const * const conf)
    :theZ(charge), theA(mass),
     theInitialZ(charge), theInitialA(mass),
     forcedTransparent(false),
     theNpInitial(0), theNnInitial(0),
     theExcitationEnergy(0.),
     initialInternalEnergy(0.),
     incomingAngularMomentum(0.,0.,0.), incomingMomentum(0.,0.,0.),
     theSpin(0.,0.,0.), theRecoilMomentum(0.,0.,0.),
     initialCenterOfMass(0.,0.,0.),
     remnant(true)
  {
    theDensity = NuclearDensityFactory::createDensity(mass, charge);

    PotentialType potentialType;
    G4bool pionPotential;
    if(conf) {
      potentialType = conf->getPotentialType();
      pionPotential = conf->getPionPotential();
    } else { // By default we don't use energy dependent
	     // potential. This is convenient for some tests.
      potentialType = IsospinPotential;
      pionPotential = true;
    }
    switch(potentialType) {
      case IsospinEnergyPotential:
        thePotential = new NuclearPotential::NuclearPotentialEnergyIsospin(theDensity, pionPotential);
        break;
      case IsospinPotential:
        thePotential = new NuclearPotential::NuclearPotentialIsospin(theDensity, pionPotential);
        break;
      case ConstantPotential:
        thePotential = new NuclearPotential::NuclearPotentialConstant(theDensity, pionPotential);
        break;
      default:
        FATAL("Unrecognized potential type at Nucleus creation." << std::endl);
        std::exit(EXIT_FAILURE);
        break;
    }

    theStore = new Store(conf);
    toBeUpdated.clear();
    blockedDelta = NULL;
  }

  Nucleus::~Nucleus() {
    delete theStore;
    delete thePotential;
    /* We don't delete the density here any more -- the Factory is caching them
    delete theDensity;*/
  }

  void Nucleus::initializeParticles()
  {
    G4INCL::ParticleType type = G4INCL::Proton;
    for(G4int i = 1; i <= theA; ++i) {
      // DEBUG("Creating particle " << i << std::endl);
      if(i == (theZ + 1)) { // Nucleons [Z+1..A] are neutrons
        type = G4INCL::Neutron;
      }

      const G4double pFermi = thePotential->getFermiMomentum(type);
      ThreeVector momentum = Random::sphereVector(pFermi);
      G4double P = momentum.mag()/pFermi;
      ThreeVector position = Random::sphereVector(theDensity->getMaxRFromP(P));
      Particle *p = new Particle(type, momentum, position);
      updatePotentialEnergy(p);
      //    p->makeParticipant(); // Force the particle to become a participant immediately (for testing)
      theStore->add(p);
    }
    initialInternalEnergy = computeTotalEnergy();
    initialCenterOfMass = computeCenterOfMass();
  }

  std::string Nucleus::dump() {
    std::stringstream ss;
    ss <<"(list ;; List of participants " << std::endl;
    ParticleList participants = theStore->getParticipants();
    for(ParticleIter i = participants.begin(); i != participants.end(); ++i) {
      ss <<"(make-particle-avatar-map " << std::endl
        << (*i)->dump() 
        << "(list ;; List of avatars in this particle" << std::endl
        << ")) ;; Close the list of avatars and the particle-avatar-map" << std::endl;
    }
    ss << ")" << std::endl;
    return ss.str();
  }

  void Nucleus::applyFinalState(FinalState *finalstate) {
    justCreated.clear();
    toBeUpdated.clear(); // Clear the list of particles to be updated by the propagation model.
    blockedDelta = NULL;
    G4double totalEnergy = 0.0;

    FinalStateValidity const validity = finalstate->getValidity();
    if(validity == ValidFS) {

      ParticleList created = finalstate->getCreatedParticles();
      for(ParticleIter iter = created.begin(); iter != created.end(); ++iter) {
        theStore->add((*iter));
        if(!(*iter)->isOutOfWell()) {
          totalEnergy += (*iter)->getEnergy() - (*iter)->getPotentialEnergy();
          justCreated.push_back((*iter));   // New particle, so we must create avatars for it
        }
      }

      ParticleList deleted = finalstate->getDestroyedParticles();
      for(ParticleIter iter = deleted.begin(); iter != deleted.end(); ++iter) {
        theStore->particleHasBeenDestroyed((*iter)->getID());
      }

      ParticleList modified = finalstate->getModifiedParticles();
      for(ParticleIter iter = modified.begin(); iter != modified.end(); ++iter) {
        theStore->particleHasBeenUpdated((*iter)->getID());
        totalEnergy += (*iter)->getEnergy() - (*iter)->getPotentialEnergy();
        toBeUpdated.push_back((*iter)); // Particle is modified so we have to create new avatars for it.
      }

      ParticleList out = finalstate->getOutgoingParticles();
      for(ParticleIter iter = out.begin(); iter != out.end(); ++iter) {
        if((*iter)->isCluster()) {
	  Cluster *clusterOut = dynamic_cast<Cluster*>((*iter));
          ParticleList const *components = clusterOut->getParticles();
          for(ParticleIter in = components->begin(); in != components->end(); ++in)
            theStore->particleHasBeenEjected((*in)->getID());
        } else {
          theStore->particleHasBeenEjected((*iter)->getID());
        }
        totalEnergy += (*iter)->getEnergy();      // No potential here because the particle is gone
        theA -= (*iter)->getA();
        theZ -= (*iter)->getZ();
        theStore->addToOutgoing(*iter);
        (*iter)->setEmissionTime(theStore->getBook()->getCurrentTime());
      }

      if(std::abs(totalEnergy - finalstate->getTotalEnergyBeforeInteraction()) > 0.1) {
        ERROR("Energy nonconservation! Energy at the beginning of the event = "
            << finalstate->getTotalEnergyBeforeInteraction()
            <<" and after interaction = "
            << totalEnergy << std::endl
            << finalstate->print());
      }
    } else if(validity == PauliBlockedFS) {
      blockedDelta = finalstate->getBlockedDelta();
    }

  }

  void Nucleus::propagateParticles(G4double /*step*/) {
    WARN("Useless Nucleus::propagateParticles -method called." << std::endl);
  }

  G4double Nucleus::computeTotalEnergy() const {
    G4double totalEnergy = 0.0;
    ParticleList inside = theStore->getParticles();
    for(ParticleIter p=inside.begin(); p!=inside.end(); ++p) {
      if((*p)->isNucleon()) // Ugly: we should calculate everything using total energies! (FIXME)
        totalEnergy += (*p)->getKineticEnergy() - (*p)->getPotentialEnergy();
      else if((*p)->isResonance()) // This is even uglier (why Proton, for instance?!)
        totalEnergy += (*p)->getEnergy() - (*p)->getPotentialEnergy() - ParticleTable::getMass(Proton);
      else
        totalEnergy += (*p)->getEnergy() - (*p)->getPotentialEnergy();
    }
    return totalEnergy;
  }

  void Nucleus::computeRecoilKinematics() {
    // If the remnant consists of only one nucleon, we need to apply a special
    // procedure to put it on mass shell.
    if(theA==1) {
      emitInsidePions();
      computeOneNucleonRecoilKinematics();
      remnant=false;
      return;
    }

    // Compute the recoil momentum and angular momentum
    theRecoilMomentum = incomingMomentum;
    theSpin = incomingAngularMomentum;

    ParticleList outgoing = theStore->getOutgoingParticles();
    for(ParticleIter p=outgoing.begin(); p!=outgoing.end(); ++p)
    {
      theRecoilMomentum -= (*p)->getMomentum();
      theSpin -= (*p)->getAngularMomentum();
    }

    // Subtract orbital angular momentum
    theCenterOfMass = computeCenterOfMass();
    theSpin -= (theCenterOfMass-initialCenterOfMass).vector(theRecoilMomentum);

    theExcitationEnergy = computeExcitationEnergy();

    G4double remnantMass = ParticleTable::getMass(theA,theZ) + theExcitationEnergy;
    theRecoilEnergy = KinematicsUtils::energy(theRecoilMomentum, remnantMass) - remnantMass;
    remnant=true;
  }

  ThreeVector Nucleus::computeCenterOfMass() const {
    ThreeVector cm(0.,0.,0.);
    G4double totalMass = 0.0;
    ParticleList inside = theStore->getParticles();
    for(ParticleIter p=inside.begin(); p!=inside.end(); ++p) {
      const G4double mass = (*p)->getMass();
      cm += (*p)->getPosition() * mass;
      totalMass += mass;
    }
    cm /= totalMass;
    return cm;
  }

  G4double Nucleus::computeExcitationEnergy() const {
    const G4double totalEnergy = computeTotalEnergy();
    const G4double separationEnergies = computeSeparationEnergyBalance();

    return totalEnergy - initialInternalEnergy - separationEnergies;
  }

  std::string Nucleus::print()
  {
    std::stringstream ss;
    ss << "Particles in the nucleus:" << std::endl
      << "Participants:" << std::endl;
    G4int counter = 1;
    ParticleList participants = theStore->getParticipants();
    for(ParticleIter p = participants.begin(); p != participants.end(); ++p) {
      ss << "index = " << counter << std::endl
        << (*p)->print();
      counter++;
    }
    ss <<"Spectators:" << std::endl;
    ParticleList spectators = theStore->getSpectators();
    for(ParticleIter p = spectators.begin(); p != spectators.end(); ++p)
      ss << (*p)->print();
    ss <<"Outgoing:" << std::endl;
    ParticleList outgoing = theStore->getOutgoingParticles();
    for(ParticleIter p = outgoing.begin(); p != outgoing.end(); ++p)
      ss << (*p)->print();

    return ss.str();
  }

  G4bool Nucleus::decayOutgoingDeltas() {
    ParticleList out = theStore->getOutgoingParticles();
    ParticleList deltas;
    for(ParticleIter i = out.begin(); i != out.end(); ++i) {
      if((*i)->isDelta()) deltas.push_back((*i));
    }
    if(deltas.empty()) return false;

    for(ParticleIter i = deltas.begin(); i != deltas.end(); ++i) {
      IAvatar *decay = new DecayAvatar((*i), 0.0, NULL);
      FinalState *fs = decay->getFinalState();
      ParticleList created = fs->getCreatedParticles();
      for(ParticleIter j = created.begin(); j != created.end(); ++j) {
          theStore->addToOutgoing(*j);
          (*j)->setEmissionTime(theStore->getBook()->getCurrentTime());
      }
      delete fs;
      delete decay;
    }
    return true;
  }

  G4bool Nucleus::decayInsideDeltas() {
    /* If there is a pion potential, do nothing (deltas will be counted as
     * excitation energy).
     * If, however, the remnant is unphysical (Z<0 or Z>A), force the deltas to
     * decay and get rid of all the pions. In case you're wondering, you can
     * end up with Z<0 or Z>A if the remnant contains more pi- than protons or
     * more pi+ than neutrons, respectively.
     */
    const G4bool unphysicalRemnant = (theZ<0 || theZ>theA);
    if(thePotential->hasPionPotential() && !unphysicalRemnant)
      return false;

    // Build a list of deltas (avoid modifying the list you are iterating on).
    ParticleList inside = theStore->getParticles();
    ParticleList deltas;
    for(ParticleIter i = inside.begin(); i != inside.end(); ++i)
      if((*i)->isDelta()) deltas.push_back((*i));

    // Loop over the deltas, make them decay
    for(ParticleIter i = deltas.begin(); i != deltas.end(); ++i) {
      // Create a forced-decay avatar. Note the last boolean parameter. Note
      // also that if the remnant is unphysical we more or less explicitly give
      // up energy conservation and CDPP by passing a NULL pointer for the
      // nucleus.
      IAvatar *decay;
      if(unphysicalRemnant)
        decay = new DecayAvatar((*i), 0.0, NULL, true);
      else
        decay = new DecayAvatar((*i), 0.0, this, true);
      FinalState *fs = decay->getFinalState();

      // The pion can be ejected only if we managed to satisfy energy
      // conservation and if pion emission does not lead to negative excitation
      // energies.
      if(fs->getValidity()==ValidFS) {
        // Apply the final state to the nucleus
        applyFinalState(fs);
      }
      delete fs;
      delete decay;
    }

    // If the remnant is unphysical, emit all the pions
    if(unphysicalRemnant) {
      DEBUG("Remnant is unphysical: Z=" << theZ << ", A=" << theA << std::endl);
      emitInsidePions();
    }

    return true;
  }

  G4bool Nucleus::decayOutgoingClusters() {
    ParticleList out = theStore->getOutgoingParticles();
    ParticleList clusters;
    for(ParticleIter i = out.begin(); i != out.end(); ++i) {
      if((*i)->isCluster()) clusters.push_back((*i));
    }
    if(clusters.empty()) return false;

    for(ParticleIter i = clusters.begin(); i != clusters.end(); ++i) {
      Cluster *cluster = dynamic_cast<Cluster*>(*i); // Can't avoid using a cast here
      ParticleList decayProducts = ClusterDecay::decay(cluster);
      for(ParticleIter j = decayProducts.begin(); j!=decayProducts.end(); ++j)
        theStore->addToOutgoing(*j);
    }
    return true;
  }

  void Nucleus::emitInsidePions() {
    DEBUG("Forcing emissions of all pions in the nucleus. This probably violates energy"
        << std::endl << "conservation (although the computation of the recoil kinematics might sweep this"
        << std::endl << "under the carpet)." << std::endl);

    // Emit the pions with this kinetic energy
    const G4double tinyPionEnergy = 0.1; // MeV

    // Push out the emitted pions
    ParticleList inside = theStore->getParticles();
    for(ParticleIter i = inside.begin(); i != inside.end(); ++i) {
      if((*i)->isPion()) {
        theZ -= (*i)->getZ();
        (*i)->setEmissionTime(theStore->getBook()->getCurrentTime());
        if((*i)->getKineticEnergy() - (*i)->getPotentialEnergy() > 0.0)
          (*i)->setEnergy((*i)->getEnergy() - (*i)->getPotentialEnergy());
        else
          (*i)->setEnergy((*i)->getMass()+tinyPionEnergy);
        (*i)->adjustMomentumFromEnergy();
        (*i)->setPotentialEnergy(0.);
        theStore->particleHasBeenEjected((*i)->getID());
        theStore->addToOutgoing(*i);
      }
    }
  }

  G4bool Nucleus::isEventTransparent() const {

    if(isForcedTransparent()) return true;

    ParticleList pL = theStore->getOutgoingParticles();
    unsigned int nIncoming = theStore->getNumberOfIncomingParticles();

    // If the number of outgoing particles is not equal to the number of
    // incoming particles, the event is not a transparent.
    if( pL.size() !=  nIncoming ) return false;

    // If any of the particles has undergone a collision, the event is not a
    // transparent.
    for(ParticleIter p = pL.begin(); p != pL.end(); ++p ) {
      if( (*p)->getNumberOfCollisions() != 0 ) return false;
      if( (*p)->getNumberOfDecays() != 0 ) return false;
    }

    return true;

  }

  void Nucleus::computeOneNucleonRecoilKinematics() {
    // We should be here only if the nucleus contains only one nucleon
// assert(theStore->getParticles().size()==1);

    DEBUG("Computing one-nucleon recoil kinematics" << std::endl);

    // No excitation energy!
    theExcitationEnergy = 0.0;

    // Move the nucleon to the outgoing list
    Particle *remN = theStore->getParticles().front();
    theA -= remN->getA();
    theZ -= remN->getZ();
    theStore->particleHasBeenEjected(remN->getID());
    theStore->addToOutgoing(remN);
    remN->setEmissionTime(theStore->getBook()->getCurrentTime());

    // Treat the special case of a remaining delta
    if(remN->isDelta()) {
      IAvatar *decay = new DecayAvatar(remN, 0.0, NULL);
      FinalState *fs = decay->getFinalState();
      // Eject the pion
      ParticleList created = fs->getCreatedParticles();
      for(ParticleIter j = created.begin(); j != created.end(); ++j)
        theStore->addToOutgoing(*j);
      delete fs;
      delete decay;
    }

    // Do different things depending on how many outgoing particles we have
    ParticleList outgoing = theStore->getOutgoingParticles();
    if(outgoing.size() == 2) {

      DEBUG("Two particles in the outgoing channel, applying exact two-body kinematics" << std::endl);

      // Can apply exact 2-body kinematics here. Keep the CM emission angle of
      // the first particle.
      Particle *p1 = outgoing.front(), *p2 = outgoing.back();
      const ThreeVector boostVector = incomingMomentum / initialEnergy;
      // Boost to the initial CM
      p1->boost(boostVector);
      const G4double sqrts = std::sqrt(initialEnergy*initialEnergy - incomingMomentum.mag2());
      const G4double pcm = KinematicsUtils::momentumInCM(sqrts, p1->getMass(), p2->getMass());
      const G4double scale = pcm/(p1->getMomentum().mag());
      // Reset the momenta
      p1->setMomentum(p1->getMomentum()*scale);
      p2->setMomentum(-p1->getMomentum());
      p1->adjustEnergyFromMomentum();
      p2->adjustEnergyFromMomentum();
      // Unboost
      p1->boost(-boostVector);
      p2->boost(-boostVector);

    } else {

      DEBUG("Trying to adjust final-state momenta to achieve energy and momentum conservation" << std::endl);

      const G4int maxIterations=8;
      G4double totalEnergy, energyScale;
      G4double val=1.E+6, oldVal=1.E+6, oldOldVal=1.E+6, oldOldOldVal;
      ThreeVector totalMomentum, deltaP;
      std::vector<ThreeVector> minMomenta;  // use it to store the particle momenta that minimize the merit function

      // Reserve the vector size
      minMomenta.reserve(outgoing.size());

      // Compute the initial total momentum
      totalMomentum.setX(0.0);
      totalMomentum.setY(0.0);
      totalMomentum.setZ(0.0);
      for(ParticleIter i=outgoing.begin(); i!=outgoing.end(); ++i)
        totalMomentum += (*i)->getMomentum();

      // Compute the initial total energy
      totalEnergy = 0.0;
      for(ParticleIter i=outgoing.begin(); i!=outgoing.end(); ++i)
        totalEnergy += (*i)->getEnergy();

      // Iterative algorithm starts here:
      for(G4int iterations=0; iterations < maxIterations; ++iterations) {

        // Save the old merit-function values
        oldOldOldVal = oldOldVal;
        oldOldVal = oldVal;
        oldVal = val;

        if(iterations%2 == 0) {
          DEBUG("Momentum step" << std::endl);
          // Momentum step: modify all the particle momenta
          deltaP = incomingMomentum - totalMomentum;
          G4double pOldTot = 0.0;
          for(ParticleIter i=outgoing.begin(); i!=outgoing.end(); ++i)
            pOldTot += (*i)->getMomentum().mag();
          for(ParticleIter i=outgoing.begin(); i!=outgoing.end(); ++i) {
            const ThreeVector mom = (*i)->getMomentum();
            (*i)->setMomentum(mom + deltaP*mom.mag()/pOldTot);
            (*i)->adjustEnergyFromMomentum();
          }
        } else {
          DEBUG("Energy step" << std::endl);
          // Energy step: modify all the particle momenta
          energyScale = initialEnergy/totalEnergy;
          for(ParticleIter i=outgoing.begin(); i!=outgoing.end(); ++i) {
            const ThreeVector mom = (*i)->getMomentum();
            G4double pScale = ((*i)->getEnergy()*energyScale - std::pow((*i)->getMass(),2))/mom.mag2();
            if(pScale>0) {
              (*i)->setEnergy((*i)->getEnergy()*energyScale);
              (*i)->adjustMomentumFromEnergy();
            }
          }
        }

        // Compute the current total momentum and energy
        totalMomentum.setX(0.0);
        totalMomentum.setY(0.0);
        totalMomentum.setZ(0.0);
        totalEnergy = 0.0;
        for(ParticleIter i=outgoing.begin(); i!=outgoing.end(); ++i) {
          totalMomentum += (*i)->getMomentum();
          totalEnergy += (*i)->getEnergy();
        }

        // Merit factor
        val = std::pow(totalEnergy - initialEnergy,2) +
          0.25*(totalMomentum - incomingMomentum).mag2();
        DEBUG("Merit function: val=" << val << ", oldVal=" << oldVal << ", oldOldVal=" << oldOldVal << ", oldOldOldVal=" << oldOldOldVal << std::endl);

        // Store the minimum
        if(val < oldVal) {
          DEBUG("New minimum found, storing the particle momenta" << std::endl);
          minMomenta.clear();
          for(ParticleIter i=outgoing.begin(); i!=outgoing.end(); ++i)
            minMomenta.push_back((*i)->getMomentum());
        }

        // Stop the algorithm if the search diverges
        if(val > oldOldVal && oldVal > oldOldOldVal) {
          DEBUG("Search is diverging, breaking out of the iteration loop: val=" << val << ", oldVal=" << oldVal << ", oldOldVal=" << oldOldVal << ", oldOldOldVal=" << oldOldOldVal << std::endl);
          break;
        }
      }

      // We should have made at least one successful iteration here
// assert(minMomenta.size()==outgoing.size());

      // Apply the optimal momenta
      DEBUG("Applying the solution" << std::endl);
      std::vector<ThreeVector>::const_iterator v = minMomenta.begin();
      for(ParticleIter i=outgoing.begin(); i!=outgoing.end(); ++i, ++v) {
        (*i)->setMomentum(*v);
        (*i)->adjustEnergyFromMomentum();
        DATABLOCK((*i)->print());
      }

    }

  }

  void Nucleus::fillEventInfo(EventInfo *eventInfo) {
    eventInfo->nParticles = 0;
    G4bool isNucleonAbsorption = false;

    G4bool isPionAbsorption = false;
    // It is possible to have pion absorption event only if the
    // projectile is pion.
    if(eventInfo->projectileType == PiPlus ||
       eventInfo->projectileType == PiMinus ||
       eventInfo->projectileType == PiZero) {
      isPionAbsorption = true;
    }

    // Outgoing particles
    ParticleList outgoingParticles = getStore()->getOutgoingParticles();

    // Check if we have a nucleon absorption event: nucleon projectile
    // and no ejected particles.
    if(outgoingParticles.size() == 0 &&
       (eventInfo->projectileType == Proton ||
	eventInfo->projectileType == Neutron)) {
      isNucleonAbsorption = true;
    }

    for( ParticleIter i = outgoingParticles.begin(); i != outgoingParticles.end(); ++i ) {
      // We have a pion absorption event only if the projectile is
      // pion and there are no ejected pions.
      if(isPionAbsorption) {
	if((*i)->isPion()) {
	  isPionAbsorption = false;
	}
      }
      eventInfo->A[eventInfo->nParticles] = (*i)->getA();
      eventInfo->Z[eventInfo->nParticles] = (*i)->getZ();
      eventInfo->emissionTime[eventInfo->nParticles] = (*i)->getEmissionTime();
      eventInfo->EKin[eventInfo->nParticles] = (*i)->getKineticEnergy();
      ThreeVector mom = (*i)->getMomentum();
      eventInfo->px[eventInfo->nParticles] = mom.getX();
      eventInfo->py[eventInfo->nParticles] = mom.getY();
      eventInfo->pz[eventInfo->nParticles] = mom.getZ();
      eventInfo->theta[eventInfo->nParticles] = Math::toDegrees(mom.theta());
      eventInfo->phi[eventInfo->nParticles] = Math::toDegrees(mom.phi());
      eventInfo->origin[eventInfo->nParticles] = -1;
      //      std::strcpy(eventInfo->history[eventInfo->nParticles],"");
      eventInfo->nParticles++;
    }
    eventInfo->nucleonAbsorption = isNucleonAbsorption;
    eventInfo->pionAbsorption = isPionAbsorption;
    eventInfo->nCascadeParticles = eventInfo->nParticles;

    // Remnant characteristics
    if(hasRemnant()) {
      eventInfo->nRemnants = 1;
      eventInfo->ARem[0] = getA();
      eventInfo->ZRem[0] = getZ();
      eventInfo->EStarRem[0] = getExcitationEnergy();
      if(eventInfo->EStarRem[0]<0.) {
	WARN("Negative excitation energy! EStarRem = " << eventInfo->EStarRem[0] << std::endl);
      }
      if(eventInfo->ARem[0]%2==0) { // even-A nucleus
	eventInfo->JRem[0] = (G4int) (getSpin().mag()/PhysicalConstants::hc + 0.5);
      } else { // odd-A nucleus
	eventInfo->JRem[0] = ((G4int) (getSpin().mag()/PhysicalConstants::hc)) + 0.5;
      }
      eventInfo->EKinRem[0] = getRecoilEnergy();
      ThreeVector mom = getRecoilMomentum();
      eventInfo->pxRem[0] = mom.getX();
      eventInfo->pyRem[0] = mom.getY();
      eventInfo->pzRem[0] = mom.getZ();
      eventInfo->thetaRem[0] = Math::toDegrees(mom.theta());
      eventInfo->phiRem[0] = Math::toDegrees(mom.phi());
    } else {
      eventInfo->nRemnants = 0;
    }

    // Global counters, flags, etc.
    eventInfo->nCollisions = getStore()->getBook()->getAcceptedCollisions();
    eventInfo->nBlockedCollisions = getStore()->getBook()->getBlockedCollisions();
    eventInfo->nDecays = getStore()->getBook()->getAcceptedDecays();
    eventInfo->nBlockedDecays = getStore()->getBook()->getBlockedDecays();
    eventInfo->firstCollisionTime = getStore()->getBook()->getFirstCollisionTime();
    eventInfo->firstCollisionXSec = getStore()->getBook()->getFirstCollisionXSec();
    eventInfo->nReflectionAvatars = getStore()->getBook()->getAvatars(SurfaceAvatarType);
    eventInfo->nCollisionAvatars = getStore()->getBook()->getAvatars(CollisionAvatarType);
    eventInfo->nDecayAvatars = getStore()->getBook()->getAvatars(DecayAvatarType);
  }

  Nucleus::ConservationBalance Nucleus::getConservationBalance(const EventInfo &theEventInfo) const {
    ConservationBalance theBalance;
    // Initialise balance variables with the incoming values
    theBalance.Z = theEventInfo.Zp + theEventInfo.Zt;
    theBalance.A = theEventInfo.Ap + theEventInfo.At;

    G4double projectileMass = 0.0;
    // FIXME: since we are not using total energies, we must set the projectile
    // mass to zero if the projectile is a pion.
    if(theEventInfo.projectileType != PiPlus &&
       theEventInfo.projectileType != PiZero &&
       theEventInfo.projectileType != PiMinus)
      projectileMass = ParticleTable::getMass(theEventInfo.projectileType);

    theBalance.energy = getInitialEnergy() - ParticleTable::getMass(theEventInfo.At, theEventInfo.Zt) - projectileMass;
    theBalance.momentum = getIncomingMomentum();

    // Process outgoing particles
    ParticleList outgoingParticles = theStore->getOutgoingParticles();
    for( ParticleIter i = outgoingParticles.begin(); i != outgoingParticles.end(); ++i ) {
      theBalance.Z -= (*i)->getZ();
      theBalance.A -= (*i)->getA();
      if((*i)->isPion()) // Ugly: we should calculate everything using total energies! (FIXME)
        theBalance.energy -= (*i)->getEnergy();
      else
        theBalance.energy -= (*i)->getKineticEnergy();
      theBalance.momentum -= (*i)->getMomentum();
    }

    theBalance.energy -= computeSeparationEnergyBalance();

    // Remnant contribution, if present
    if(hasRemnant()) {
      theBalance.Z -= getZ();
      theBalance.A -= getA();
      theBalance.energy -= //ParticleTable::getMass(getA(),getZ()) +
        getExcitationEnergy() + getRecoilEnergy();
      theBalance.momentum -= getRecoilMomentum();
    }

    return theBalance;
  }

}

#endif
