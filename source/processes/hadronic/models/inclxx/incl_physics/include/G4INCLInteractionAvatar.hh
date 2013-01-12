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

/* \file G4INCLInteractionAvatar.hh
 * \brief Virtual class for interaction avatars.
 *
 * This class is inherited by decay and collision avatars. The goal is to
 * provide a uniform treatment of common physics, such as Pauli blocking,
 * enforcement of energy conservation, etc.
 *
 *  \date Mar 1st, 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLINTERACTIONAVATAR_HH_
#define G4INCLINTERACTIONAVATAR_HH_

#include "G4INCLIAvatar.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLRootFinder.hh"
#include "G4INCLKinematicsUtils.hh"

namespace G4INCL {

  class InteractionAvatar : public G4INCL::IAvatar {
    public:
      InteractionAvatar(G4double, G4INCL::Nucleus*, G4INCL::Particle*);
      InteractionAvatar(G4double, G4INCL::Nucleus*, G4INCL::Particle*, G4INCL::Particle*);
      virtual ~InteractionAvatar();

      /// \brief Target accuracy in the determination of the local-energy Q-value
      static const G4double locEAccuracy;
      /// \brief Max number of iterations for the determination of the local-energy Q-value
      static const G4int maxIterLocE;

    protected:
      virtual G4INCL::IChannel* getChannel() const = 0;

      G4bool bringParticleInside(Particle * const p);

      /** \brief Apply local-energy transformation, if appropriate
       *
       * \param p particle to apply the transformation to
       */
      void preInteractionLocalEnergy(Particle * const p);

      /** \brief Store the state of the particles before the interaction
       *
       * If the interaction cannot be realised for any reason, we will need to
       * restore the particle state as it was before. This is done by calling
       * the restoreParticles() method.
       */
      void preInteractionBlocking();

      void preInteraction();
      FinalState *postInteraction(FinalState *);

      /** \brief Restore the state of both particles.
       *
       * The state must first be stored by calling preInteractionBlocking().
       */
      void restoreParticles() const;

      /// \brief true if the given avatar should use local energy
      G4bool shouldUseLocalEnergy() const {
        if(!theNucleus) return false;
        LocalEnergyType theLocalEnergyType;
        if(getType()==DecayAvatarType || isPiN)
          theLocalEnergyType = theNucleus->getStore()->getConfig()->getLocalEnergyPiType();
        else
          theLocalEnergyType = theNucleus->getStore()->getConfig()->getLocalEnergyBBType();

        const G4bool firstAvatar = (theNucleus->getStore()->getBook()->getAcceptedCollisions() == 0);
        return ((theLocalEnergyType == FirstCollisionLocalEnergy && firstAvatar) ||
            theLocalEnergyType == AlwaysLocalEnergy);
      }

      G4INCL::Nucleus *theNucleus;
      G4INCL::Particle *particle1, *particle2;
      ThreeVector boostVector;
      ParticleType oldParticle1Type, oldParticle2Type;
      G4double oldParticle1Energy, oldParticle2Energy, oldTotalEnergy, oldXSec;
      G4double oldParticle1Potential, oldParticle2Potential;
      G4double oldParticle1Mass, oldParticle2Mass;
      G4double oldParticle1Helicity, oldParticle2Helicity;
      ThreeVector oldParticle1Momentum, oldParticle2Momentum;
      ThreeVector oldParticle1Position, oldParticle2Position;
      G4bool isPiN;

    private:
      /// \brief RootFunctor-derived object for enforcing energy conservation in N-N.
      class ViolationEMomentumFunctor : public RootFunctor {
        public:
          /** \brief Prepare for calling the () operator and scaleParticleMomenta
           *
           * The constructor sets the private class members.
           */
          ViolationEMomentumFunctor(Nucleus * const nucleus, FinalState const * const finalState, ThreeVector const * const boost, const G4bool localE);
          virtual ~ViolationEMomentumFunctor() { particleMomenta.clear(); }

          /** \brief Compute the energy-conservation violation.
           *
           * \param x scale factor for the particle momenta
           * \return the energy-conservation violation
           */
          G4double operator()(const G4double x) const;

          /// \brief Clean up after root finding
          void cleanUp(const G4bool success) const;

        private:
          /// \brief List of final-state particles.
          ParticleList finalParticles;
          /// \brief CM particle momenta, as determined by the channel.
          std::list<ThreeVector> particleMomenta;
          /// \brief Total energy before the interaction.
          G4double initialEnergy;
          /// \brief Pointer to the nucleus
          Nucleus *theNucleus;
          /// \brief Pointer to the boost vector
          ThreeVector const *boostVector;
          /// \brief true if we must apply local energy to nucleons
          G4bool hasLocalEnergy;
          /// \brief true if we must apply local energy to deltas
          G4bool hasLocalEnergyDelta;

          /// \brief True if we should use local energy
          const G4bool shouldUseLocalEnergy;

          /** \brief Scale the momenta of the modified and created particles.
           *
           * Set the momenta of the modified and created particles to alpha times
           * their original momenta (stored in particleMomenta). You must call
           * init() before using this method.
           * 
           * \param alpha scale factor
           */
          void scaleParticleMomenta(const G4double alpha) const;

      };

      /// \brief RootFunctor-derived object for enforcing energy conservation in pi-N.
      class ViolationEEnergyFunctor : public RootFunctor {
        public:
          /** \brief Prepare for calling the () operator and scaleParticleMomenta
           *
           * The constructor sets the private class members.
           */
          ViolationEEnergyFunctor(Nucleus * const nucleus, FinalState const * const finalState);
          virtual ~ViolationEEnergyFunctor() {}

          /** \brief Compute the energy-conservation violation.
           *
           * \param x scale factor for the particle momenta
           * \return the energy-conservation violation
           */
          G4double operator()(const G4double x) const;

          /// \brief Clean up after root finding
          void cleanUp(const G4bool success) const;

          /** \brief Set the energy of the particle.
           *
           * \param energy
           */
          void setParticleEnergy(const G4double energy) const;

        private:
          /// \brief Total energy before the interaction.
          G4double initialEnergy;
          /// \brief Pointer to the nucleus.
          Nucleus *theNucleus;
          /// \brief The final-state particle.
          Particle *theParticle;
          /// \brief The initial energy of the particle.
          G4double theEnergy;
          /// \brief The initial momentum of the particle.
          ThreeVector theMomentum;
          /** \brief Threshold for the energy of the particle
           *
           * The particle (a delta) cannot have less than this energy.
           */
          G4double energyThreshold;
      };

      RootFunctor *violationEFunctor;

    protected:
      /** \brief Enforce energy conservation.
       *
       * Final states generated by the channels might violate energy conservation
       * because of different reasons (energy-dependent potentials, local
       * energy...). This conservation law must therefore be enforced by hand. We
       * do so by rescaling the momenta of the final-state particles in the CM
       * frame. If this turns out to be impossible, this method returns false.
       *
       * \return true if the algorithm succeeded
       */
      G4bool enforceEnergyConservation(FinalState * const fs);

  };

}

#endif /* G4INCLINTERACTIONAVATAR_HH_ */
