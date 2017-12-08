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

#ifndef G4INCLCascade_hh
#define G4INCLCascade_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLCascadeAction.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLGlobalInfo.hh"
#include "G4INCLLogger.hh"
#include "G4INCLConfig.hh"
#include "G4INCLRootFinder.hh"

namespace G4INCL {
  class INCL {
    public:
      INCL(Config const * const config);

      ~INCL();

      /// \brief Dummy copy constructor to silence Coverity warning
      INCL(const INCL &rhs);

      /// \brief Dummy assignment operator to silence Coverity warning
      INCL &operator=(const INCL &rhs);

      G4bool prepareReaction(const ParticleSpecies &projectileSpecies, const G4double kineticEnergy, const G4int A, const G4int Z, const G4int S);
      G4bool initializeTarget(const G4int A, const G4int Z, const G4int S);
      inline const EventInfo &processEvent() {
        return processEvent(
            theConfig->getProjectileSpecies(),
            theConfig->getProjectileKineticEnergy(),
            theConfig->getTargetA(),
            theConfig->getTargetZ(),
            theConfig->getTargetS()
            );
      }
      const EventInfo &processEvent(
          ParticleSpecies const &projectileSpecies,
          const G4double kineticEnergy,
          const G4int targetA,
          const G4int targetZ,
          const G4int targetS
          );

      void finalizeGlobalInfo(Random::SeedVector const &initialSeeds);
      const GlobalInfo &getGlobalInfo() const { return theGlobalInfo; }
      

    private:
      IPropagationModel *propagationModel;
      G4int theA, theZ, theS;
      G4bool targetInitSuccess;
      G4double maxImpactParameter;
      G4double maxUniverseRadius;
      G4double maxInteractionDistance;
      G4double fixedImpactParameter;
      CascadeAction *cascadeAction;
      Config const * const theConfig;
      Nucleus *nucleus;
      G4bool forceTransparent;
      
      EventInfo theEventInfo;
      GlobalInfo theGlobalInfo;

      /// \brief Remnant size below which cascade stops
      G4int minRemnantSize;

      /// \brief Class to adjust remnant recoil
      class RecoilFunctor : public RootFunctor {
        public:
          /** \brief Prepare for calling the () operator and scaleParticleEnergies
           *
           * The constructor sets the private class members.
           */
          RecoilFunctor(Nucleus * const n, const EventInfo &ei) :
            RootFunctor(0., 1E6),
            nucleus(n),
            outgoingParticles(n->getStore()->getOutgoingParticles()),
            theEventInfo(ei) {
              for(ParticleIter p=outgoingParticles.begin(), e=outgoingParticles.end(); p!=e; ++p) {
                particleMomenta.push_back((*p)->getMomentum());
                particleKineticEnergies.push_back((*p)->getKineticEnergy());
              }
              ProjectileRemnant * const aPR = n->getProjectileRemnant();
              if(aPR && aPR->getA()>0) {
                particleMomenta.push_back(aPR->getMomentum());
                particleKineticEnergies.push_back(aPR->getKineticEnergy());
                outgoingParticles.push_back(aPR);
              }
            }
          virtual ~RecoilFunctor() {}

          /** \brief Compute the energy-conservation violation.
           *
           * \param x scale factor for the particle energies
           * \return the energy-conservation violation
           */
          G4double operator()(const G4double x) const {
            scaleParticleEnergies(x);
            return nucleus->getConservationBalance(theEventInfo,true).energy;
          }

          /// \brief Clean up after root finding
          void cleanUp(const G4bool success) const {
            if(!success)
              scaleParticleEnergies(1.);
          }

        private:
          /// \brief Pointer to the nucleus
          Nucleus *nucleus;
          /// \brief List of final-state particles.
          ParticleList outgoingParticles;
          // \brief Reference to the EventInfo object
          EventInfo const &theEventInfo;
          /// \brief Initial momenta of the outgoing particles
          std::list<ThreeVector> particleMomenta;
          /// \brief Initial kinetic energies of the outgoing particles
          std::list<G4double> particleKineticEnergies;

          /** \brief Scale the kinetic energies of the outgoing particles.
           *
           * \param rescale scale factor
           */
          void scaleParticleEnergies(const G4double rescale) const {
            // Rescale the energies (and the momenta) of the outgoing particles.
            ThreeVector pBalance = nucleus->getIncomingMomentum();
            std::list<ThreeVector>::const_iterator iP = particleMomenta.begin();
            std::list<G4double>::const_iterator iE = particleKineticEnergies.begin();
            for( ParticleIter i = outgoingParticles.begin(), e = outgoingParticles.end(); i!=e; ++i, ++iP, ++iE)
            {
              const G4double mass = (*i)->getMass();
              const G4double newKineticEnergy = (*iE) * rescale;

              (*i)->setMomentum(*iP);
              (*i)->setEnergy(mass + newKineticEnergy);
              (*i)->adjustMomentumFromEnergy();

              pBalance -= (*i)->getMomentum();
            }

            nucleus->setMomentum(pBalance);
            const G4double remnantMass = ParticleTable::getTableMass(nucleus->getA(),nucleus->getZ()) + nucleus->getExcitationEnergy();
            const G4double pRem2 = pBalance.mag2();
            const G4double recoilEnergy = pRem2/
              (std::sqrt(pRem2+remnantMass*remnantMass) + remnantMass);
            nucleus->setEnergy(remnantMass + recoilEnergy);
          }
      };

      /// \brief Class to adjust remnant recoil in the reaction CM system
      class RecoilCMFunctor : public RootFunctor {
        public:
          /** \brief Prepare for calling the () operator and scaleParticleEnergies
           *
           * The constructor sets the private class members.
           */
          RecoilCMFunctor(Nucleus * const n, const EventInfo &ei) :
            RootFunctor(0., 1E6),
            nucleus(n),
            theIncomingMomentum(nucleus->getIncomingMomentum()),
            outgoingParticles(n->getStore()->getOutgoingParticles()),
            theEventInfo(ei) {
              thePTBoostVector = nucleus->getIncomingMomentum()/nucleus->getInitialEnergy();
              for(ParticleIter p=outgoingParticles.begin(), e=outgoingParticles.end(); p!=e; ++p) {
                (*p)->boost(thePTBoostVector);
                particleCMMomenta.push_back((*p)->getMomentum());
              }
              ProjectileRemnant * const aPR = n->getProjectileRemnant();
              if(aPR && aPR->getA()>0) {
                aPR->boost(thePTBoostVector);
                particleCMMomenta.push_back(aPR->getMomentum());
                outgoingParticles.push_back(aPR);
              }
            }
          virtual ~RecoilCMFunctor() {}

          /** \brief Compute the energy-conservation violation.
           *
           * \param x scale factor for the particle energies
           * \return the energy-conservation violation
           */
          G4double operator()(const G4double x) const {
            scaleParticleCMMomenta(x);
            return nucleus->getConservationBalance(theEventInfo,true).energy;
          }

          /// \brief Clean up after root finding
          void cleanUp(const G4bool success) const {
            if(!success)
              scaleParticleCMMomenta(1.);
          }

        private:
          /// \brief Pointer to the nucleus
          Nucleus *nucleus;
          /// \brief Projectile-target CM boost vector
          ThreeVector thePTBoostVector;
          /// \brief Incoming momentum
          ThreeVector theIncomingMomentum;
          /// \brief List of final-state particles.
          ParticleList outgoingParticles;
          // \brief Reference to the EventInfo object
          EventInfo const &theEventInfo;
          /// \brief Initial CM momenta of the outgoing particles
          std::list<ThreeVector> particleCMMomenta;

          /** \brief Scale the kinetic energies of the outgoing particles.
           *
           * \param rescale scale factor
           */
          void scaleParticleCMMomenta(const G4double rescale) const {
            // Rescale the CM momenta of the outgoing particles.
            ThreeVector remnantMomentum = theIncomingMomentum;
            std::list<ThreeVector>::const_iterator iP = particleCMMomenta.begin();
            for( ParticleIter i = outgoingParticles.begin(), e = outgoingParticles.end(); i!=e; ++i, ++iP)
            {
              (*i)->setMomentum(*iP * rescale);
              (*i)->adjustEnergyFromMomentum();
              (*i)->boost(-thePTBoostVector);

              remnantMomentum -= (*i)->getMomentum();
            }

            nucleus->setMomentum(remnantMomentum);
            const G4double remnantMass = ParticleTable::getTableMass(nucleus->getA(),nucleus->getZ()) + nucleus->getExcitationEnergy();
            const G4double pRem2 = remnantMomentum.mag2();
            const G4double recoilEnergy = pRem2/
              (std::sqrt(pRem2+remnantMass*remnantMass) + remnantMass);
            nucleus->setEnergy(remnantMass + recoilEnergy);
          }
      };

      /** \brief Rescale the energies of the outgoing particles.
       *
       * Allow for the remnant recoil energy by rescaling the energy (and
       * momenta) of the outgoing particles.
       */
      void rescaleOutgoingForRecoil();

#ifndef INCLXX_IN_GEANT4_MODE
      /** \brief Run global conservation checks
       *
       * Check that energy and momentum are correctly conserved. If not, issue
       * a warning.
       *
       * Also feeds the balance variables in theEventInfo.
       *
       * \param afterRecoil whether to take into account nuclear recoil
       */
      void globalConservationChecks(G4bool afterRecoil);
#endif

      /** \brief Stopping criterion for the cascade
       *
       * Returns true if the cascade should continue, and false if any of the
       * stopping criteria is satisfied.
       */
      G4bool continueCascade();

      /** \brief Make a projectile pre-fragment out of geometrical spectators
       *
       * The projectile pre-fragment is assigned an excitation energy given
       * by \f$E_\mathrm{sp}-E_\mathrm{i,A}\f$, where \f$E_\mathrm{sp}\f$ is the
       * sum of the energies of the spectator particles, and \f$E_\mathrm{i,A}\f$
       * is the sum of the smallest \f$A\f$ particle energies initially present
       * in the projectile, \f$A\f$ being the mass of the projectile
       * pre-fragment. This is equivalent to assuming that the excitation
       * energy is given by the sum of the transitions of all excited
       * projectile components to the "holes" left by the participants.
       *
       * This method can modify the outgoing list and adds a projectile
       * pre-fragment.
       *
       * \return the number of dynamical spectators that were merged back in
       *         the projectile
       */
      G4int makeProjectileRemnant();

      /** \brief Make a compound nucleus
       *
       * Selects the projectile components that can actually enter their
       * potential and puts them into the target nucleus. If the CN excitation
       * energy turns out to be negative, the event is considered a
       * transparent. This method modifies theEventInfo and theGlobalInfo.
       */
      void makeCompoundNucleus();

      /// \brief Initialise the cascade
      G4bool preCascade(ParticleSpecies const &projectileSpecies, const G4double kineticEnergy);

      /// \brief The actual cascade loop
      void cascade();

      /// \brief Finalise the cascade and clean up
      void postCascade();

      /** \brief Initialise the maximum interaction distance.
       *
       * Used in forced CN events.
       */
      void initMaxInteractionDistance(ParticleSpecies const &p, const G4double kineticEnergy);

      /** \brief Initialize the universe radius.
       *
       * Used for determining the energy-dependent size of the volume particles
       * live in.
       */
      void initUniverseRadius(ParticleSpecies const &p, const G4double kineticEnergy, const G4int A, const G4int Z);

      /// \brief Update global counters and other members of theGlobalInfo object
      void updateGlobalInfo();
  };
}

#endif
