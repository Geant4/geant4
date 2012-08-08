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
// INCL++ revision: v5.1.2
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4INCLCascade_hh
#define G4INCLCascade_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLEventAction.hh"
#include "G4INCLPropagationAction.hh"
#include "G4INCLAvatarAction.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLGlobalInfo.hh"
#include "G4INCLLogger.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {
  class INCL {
    public:
      INCL(Config const * const config);

      ~INCL();

      G4bool prepareReaction(const ParticleSpecies &projectileSpecies, const G4double kineticEnergy, const G4int A, const G4int Z);
      G4bool initializeTarget(const G4int A, const G4int Z);
      inline const EventInfo &processEvent() {
        return processEvent(
            theConfig->getProjectileSpecies(),
            theConfig->getProjectileKineticEnergy(),
            theConfig->getTargetA(),
            theConfig->getTargetZ()
            );
      }
      const EventInfo &processEvent(
          ParticleSpecies const &projectileSpecies,
          const G4double kineticEnergy,
          const G4int targetA,
          const G4int targetZ
          );

      void finaliseGlobalInfo();
      const GlobalInfo &getGlobalInfo() const { return theGlobalInfo; }

      std::string configToString() { return theConfig->echo(); }

    private:
      IPropagationModel *propagationModel;
      G4int theA, theZ;
      G4bool targetInitSuccess;
      G4double maxImpactParameter;
      G4double maxUniverseRadius;
      G4double maxInteractionDistance;
      G4double fixedImpactParameter;
      EventAction *eventAction;
      PropagationAction *propagationAction;
      AvatarAction *avatarAction;
      Config const * const theConfig;
      Nucleus *nucleus;

      EventInfo theEventInfo;
      GlobalInfo theGlobalInfo;

      /// \brief Remnant size below which cascade stops
      static const G4int minRemnantSize;

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
       * \param afterRescaling whether to take into account nuclear recoil
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
       * The projectile pre-fragment is assigned an excitation energy given by
       * \f$E_\mathrm{sp}-E_\mathrm{i,A}\f$, where \f$E_\mathrm{sp}\f$ is the
       * sum of the energies of the spectator particles, and
       * \f$E_\mathrm{i,A}\f$ is the sum of the smallest \f$A\f$ particle
       * energies initially present in the projectile, \f$A\f$ being the mass
       * of the projectile pre-fragment. This is equivalent to assuming that
       * the excitation energy is given by the sum of the transitions of all
       * excited projectile components to the "holes" left by the participants.
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
      G4bool preCascade(ParticleSpecies const projectileSpecies, const G4double kineticEnergy);

      /// \brief The actual cascade loop
      void cascade();

      /// \brief Finalise the cascade and clean up
      void postCascade();

      /** \brief Initialise the maximum interaction distance.
       *
       * Used in forced CN events.
       */
      void initMaxInteractionDistance(Cluster const * const c);

      /** \brief Initialize the universe radius.
       *
       * Used for determining the energy-dependent size of the volume particles
       * live in.
       */
      void initUniverseRadius(ParticleSpecies const &p, const G4double kineticEnergy, const G4int A, const G4int Z);
  };
}

#endif
