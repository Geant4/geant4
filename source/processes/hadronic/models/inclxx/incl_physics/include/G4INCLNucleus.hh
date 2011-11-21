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
 * G4INCLNucleus.hh
 *
 *  Created on: Jun 5, 2009
 *      Author: Pekka Kaitaniemi
 */

#ifndef G4INCLNUCLEUS_HH_
#define G4INCLNUCLEUS_HH_

#include <list>
#include <string>

#include "G4INCLParticle.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLCluster.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLStore.hh"
#include "G4INCLNuclearDensity.hh"
#include "G4INCLINuclearPotential.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLConfig.hh"
#include "G4INCLConfigEnums.hh"

namespace G4INCL {

  class Nucleus {
  public:
    Nucleus(G4int mass, G4int charge, Config const * const conf);
    virtual ~Nucleus();

    /**
     * Generate the initial distribution of particles. At the beginning
     * all particles are assigned as spectators.
     */
    void initializeParticles();

    /**
     * Insert a new participant (e.g. a projectile) to the nucleus.
     */
    void insertParticipant(Particle *p) {
      p->makeParticipant(); // The projectile particle is a participant
      theZ += p->getZ();
      theA += p->getA();
      theStore->particleHasEntered(p);
      if(p->isNucleon()) {
        theNpInitial += Math::heaviside(ParticleTable::getIsospin(p->getType()));
        theNnInitial += Math::heaviside(-ParticleTable::getIsospin(p->getType()));
      }
    };

    /**
     * Calculate the transmission probability for particle p
     */
    G4double getTransmissionProbability(Particle const * const p);

    /**
     * Apply reaction final state information to the nucleus.
     */
    void applyFinalState(FinalState *);

    G4int getA() const { return theA; };
    G4int getZ() const { return theZ; };

    G4int getInitialA() const { return theInitialA; };
    G4int getInitialZ() const { return theInitialZ; };

    /**
     * Get the list of particles that were created by the last applied final state
     */
    ParticleList const &getCreatedParticles() const { return justCreated; }

    /**
     * Get the list of particles that were updated by the last applied final state
     */
    ParticleList const &getUpdatedParticles() const { return toBeUpdated; }

    /// \brief Get the delta that could not decay
    Particle *getBlockedDelta() const { return blockedDelta; }

    /**
     * Propagate the particles one time step.
     *
     * @param step length of the time step
     */
    void propagateParticles(G4double step);

    G4int getNumberOfProjectileProtons() const { return theNpInitial; };
    G4int getNumberOfProjectileNeutrons() const { return theNnInitial; };

    /** \brief Outgoing - incoming separation energies.
     *
     * Used by CDPP.
     */
    G4double computeSeparationEnergyBalance() const {
      G4double S = 0.0;
      ParticleList outgoing = theStore->getOutgoingParticles();
      for(ParticleIter i = outgoing.begin(); i != outgoing.end(); ++i)
        if((*i)->isNucleon() || (*i)->isResonance())
          S += ParticleTable::getSeparationEnergy((*i)->getType());
        else if((*i)->isCluster()) {
          S += (*i)->getZ() * ParticleTable::getSeparationEnergy(Proton)
            + ((*i)->getA() - (*i)->getZ()) * ParticleTable::getSeparationEnergy(Neutron);
        }

      S -= theNpInitial * ParticleTable::getSeparationEnergy(Proton);
      S -= theNnInitial * ParticleTable::getSeparationEnergy(Neutron);
      return S;
    }

    /** \brief Force the decay of outgoing deltas.
     *
     * \return true if any delta was forced to decay.
     */
    G4bool decayOutgoingDeltas();

    /** \brief Force the decay of deltas inside the nucleus.
     *
     * \return true if any delta was forced to decay.
     */
    G4bool decayInsideDeltas();

    /** \brief Force the decay of unstable outgoing clusters.
     *
     * \return true if any cluster was forced to decay.
     */
    G4bool decayOutgoingClusters();

    /// \brief Force emission of all pions inside the nucleus.
    void emitInsidePions();

    /** \brief Compute the recoil momentum and spin of the nucleus. */
    void computeRecoilKinematics();

    /** \brief Compute the current center-of-mass position.
     *
     * \return the center-of-mass position vector [fm].
     */
    ThreeVector computeCenterOfMass() const;

    /** \brief Compute the current total energy.
     *
     * \return the total energy [MeV]
     */
    G4double computeTotalEnergy() const;

    /** \brief Compute the current excitation energy.
     *
     * \return the excitation energy [MeV]
     */
    G4double computeExcitationEnergy() const;

    /** \brief Set the incoming angular-momentum vector. */
    void setIncomingAngularMomentum(const ThreeVector &j) {
      incomingAngularMomentum = j;
    }

    /** \brief Set the incoming momentum vector. */
    void setIncomingMomentum(const ThreeVector &p) {
      incomingMomentum = p;
    }

    /** \brief Get the incoming momentum vector. */
    const ThreeVector &getIncomingMomentum() const {
      return incomingMomentum;
    }

    /** \brief Set the initial energy. */
    void setInitialEnergy(const G4double e) { initialEnergy = e; }

    /** \brief Get the initial energy. */
    G4double getInitialEnergy() const { return initialEnergy; }

    /** \brief Get the recoil energy of the nucleus.
     *
     * Method computeRecoilKinematics() should be called first.
     */
    G4double getRecoilEnergy() const { return theRecoilEnergy; }

    /** \brief Get the excitation energy of the nucleus.
     *
     * Method computeRecoilKinematics() should be called first.
     */
    G4double getExcitationEnergy() const { return theExcitationEnergy; }

    /** \brief Get the spin of the nucleus.
     *
     * Method computeRecoilKinematics() should be called first.
     */
    ThreeVector const &getSpin() const { return theSpin; }

    /** \brief Get the recoil momentum of the nucleus.
     *
     * Method computeRecoilKinematics() should be called first.
     */
    ThreeVector const &getRecoilMomentum() const { return theRecoilMomentum; }

    /** \brief Set the recoil momentum of the nucleus
     *
     * Can be used to override the recoil momentum computed by
     * computeRecoilKinematics();
     * */
    void setRecoilMomentum(const ThreeVector &p) { theRecoilMomentum = p; }

    /** \brief Set the recoil energy of the nucleus
     *
     * Can be used to override the recoil energy computed by
     * computeRecoilKinematics();
     * */
    void setRecoilEnergy(G4double energy) { theRecoilEnergy = energy; }

    /**
     * Mark a particle as a participant.
     *
     * @param p poG4inter to a particle
     */
    void participate(G4INCL::Particle *p);

    NuclearDensity* getDensity() const { return theDensity; };
    NuclearPotential::INuclearPotential* getPotential() const { return thePotential; };

    /// \brief Update the particle potential energy.
    inline void updatePotentialEnergy(G4INCL::Particle *p) {
      p->setPotentialEnergy(thePotential->computePotentialEnergy(p));
    }

    ///\brief Returns true if the nucleus contains any deltas.
    inline G4bool containsDeltas() {
      ParticleList inside = theStore->getParticles();
      for(ParticleIter i=inside.begin(); i!=inside.end(); ++i)
        if((*i)->isDelta()) return true;
      return false;
    }

    /** \brief Modify particle that enters the nucleus.
     *
     * Modify the particle momentum and/or position when the particle enters
     * the nucleus.
     *
     * \param particle poG4inter to entering particle
     * \return poG4inter to modified particle
     */
    Particle *particleEnters(Particle *particle);

    /** \brief Modify particle that leaves the nucleus.
     *
     * Modify the particle momentum and/or position when the particle leaves
     * the nucleus.
     *
     * \param particle poG4inter to leaving particle
     * \return poG4inter to modified particle
     */
    Particle *particleLeaves(Particle *particle);

    /** \brief Get the maximum allowed radius for a given particle.
     * 
     * Calls the NuclearDensity::getMaxRFromP() method for nucleons and deltas,
     * and the NuclearDensity::getTrasmissionRadius() method for pions.
     *
     * \param particle poG4inter to a particle
     * \return surface radius
     */
    G4double getSurfaceRadius(Particle const * const particle) const {
      if(particle->isPion())
        // Temporarily set RPION = RMAX
        return theDensity->getMaximumRadius();
        //return 0.5*(theDensity->getTransmissionRadius(particle)+theDensity->getMaximumRadius());
      else {
        const G4double pr = particle->getMomentum().mag()/thePotential->getFermiMomentum(particle);
        return theDensity->getMaxRFromP(pr);
      }
    }

    /**
     * PrG4int the nucleus info
     */
    std::string prG4int();

    std::string dump();

    Store* getStore() const {return theStore; };
    void setStore(Store *s) {
      delete theStore;
      theStore = s;
    };

    G4double getInitialInternalEnergy() const { return initialInternalEnergy; };

    /** \brief Is the event transparent?
     *
     * To be called at the end of the cascade.
     **/
    G4bool isEventTransparent() const;

    /** \brief Does the nucleus give a cascade remnant?
     *
     * To be called after computeRecoilKinematics().
     **/
    G4bool hasRemnant() const { return remnant; }

    void forceTransparent() { forcedTransparent=true; }
    G4bool isForcedTransparent() const { return forcedTransparent; }

    /**
     * Fill the event info which contains INCL output data
     */
    //    void fillEventInfo(Results::EventInfo *eventInfo);
    void fillEventInfo(EventInfo *eventInfo);

  private:
    /** \brief Compute the recoil kinematics for a 1-nucleon remnant.
     *
     * Puts the remnant nucleon on mass shell and tries to enforce approximate
     * energy conservation by modifying the masses of the outgoing particles.
     */
    void computeOneNucleonRecoilKinematics();

  private:
    G4int theZ, theA;
    G4int theInitialZ, theInitialA;
    G4bool forcedTransparent;
    G4int theNpInitial, theNnInitial;
    G4double theExcitationEnergy;
    G4double initialInternalEnergy;
    ThreeVector incomingAngularMomentum, incomingMomentum;
    ThreeVector theSpin, theRecoilMomentum, theCenterOfMass;
    ThreeVector initialCenterOfMass;
    G4bool remnant;

    ParticleList toBeUpdated;
    ParticleList justCreated;
    Particle *blockedDelta;
    NuclearDensity *theDensity;
    NuclearPotential::INuclearPotential *thePotential;
    G4double theRecoilEnergy;
    G4double initialEnergy;
    Store *theStore;
  };

}

#endif /* G4INCLNUCLEUS_HH_ */
