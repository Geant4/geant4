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

#include "G4INCLDecayAvatar.hh"

#include "G4INCLDeltaDecayChannel.hh"
#include "G4INCLPauliBlocking.hh"
#include <sstream>
#include <string>
// #include <cassert>

namespace G4INCL {

  DecayAvatar::DecayAvatar(G4INCL::Particle *aParticle, G4double time, G4INCL::Nucleus *n, G4bool force)
    : InteractionAvatar(time, n, aParticle), forced(force),
      incidentDirection(aParticle->getMomentum())
  {
    setType(DecayAvatarType);
  }

  DecayAvatar::~DecayAvatar() {

  }

  G4INCL::IChannel* DecayAvatar::getChannel() const
  {
    if(particle1->isDelta()) {
      DEBUG("DeltaDecayChannel chosen." << std::endl);
      return new DeltaDecayChannel(theNucleus, particle1, incidentDirection);
    }
    else
      return NULL;
  }

  void DecayAvatar::preInteraction() {
    InteractionAvatar::preInteraction();
  }

  FinalState *DecayAvatar::postInteraction(FinalState *fs) {
    // Make sure we have at least two particles in the final state
// assert(fs->getModifiedParticles().size() + fs->getCreatedParticles().size() - fs->getDestroyedParticles().size() >= 2);

    if(!forced) { // Normal decay

      // Call the postInteraction method of the parent class
      // (provides Pauli blocking and enforces energy conservation)
      fs = InteractionAvatar::postInteraction(fs);

      if(fs->getValidity() == PauliBlockedFS)
        /* If the decay was Pauli-blocked, make sure the propagation model
         * generates a new decay avatar on the next call to propagate().
         *
         * \bug{Note that we don't generate new decay avatars for deltas that
         * could not satisfy energy conservation. This is in keeping with
         * INCL4.6, but doesn't seem to make much sense to me (DM), as energy
         * conservation can be impossible to satisfy due to weird local-energy
         * conditions, for example, that evolve with time.}
         */
        fs->setBlockedDelta(particle1);

    } else { // Forced decay
      ParticleList created = fs->getCreatedParticles();

      // Try to enforce energy conservation
      fs->setTotalEnergyBeforeInteraction(oldTotalEnergy);
      const G4bool success = enforceEnergyConservation(fs);
      if(!success) {
        DEBUG("Enforcing energy conservation: failed!" << std::endl);

        if(theNucleus) {
          // Restore the state of the initial particles
          restoreParticles();

          // Delete newly created particles
          for( ParticleIter i = created.begin(); i != created.end(); ++i )
            delete *i;

          FinalState *fsBlocked = new FinalState;
          delete fs;
          fsBlocked->makeNoEnergyConservation();
          fsBlocked->setTotalEnergyBeforeInteraction(0.0);

          return fsBlocked; // Interaction is blocked. Return an empty final state.
        } else {
          // If there is no nucleus we have to continue anyway, even if energy
          // conservation failed. We cannot afford producing unphysical
          // remnants.
          DEBUG("No nucleus, continuing anyway." << std::endl);
        }
      } else {
        DEBUG("Enforcing energy conservation: success!" << std::endl);
      }

      if(theNucleus) {
        ParticleList modified = fs->getModifiedParticles();

        // Copy the final state, but don't include the pion (as if it had been
        // emitted right away).
        FinalState *emissionFS = new FinalState;
        for(ParticleIter i=modified.begin(); i!=modified.end(); ++i)
          emissionFS->addModifiedParticle(*i);

        // Test CDPP blocking
        G4bool isCDPPBlocked = Pauli::isCDPPBlocked(created, theNucleus);

        if(isCDPPBlocked) {
          DEBUG("CDPP: Blocked!" << std::endl);

          // Restore the state of both particles
          restoreParticles();

          // Delete newly created particles
          for( ParticleIter i = created.begin(); i != created.end(); ++i )
            delete *i;

          FinalState *fsBlocked = new FinalState;
          delete fs;
	  delete emissionFS;

          fsBlocked->makePauliBlocked();
          fsBlocked->setTotalEnergyBeforeInteraction(0.0);

          return fsBlocked; // Interaction is blocked. Return an empty final state.
        }
        DEBUG("CDPP: Allowed!" << std::endl);

        // If all went well (energy conservation enforced and CDPP satisfied),
        // delete the auxiliary final state
        delete emissionFS;

      }
    }

    // If there is a nucleus, increment the counters
    if(theNucleus) {
      switch(fs->getValidity()) {
        case PauliBlockedFS:
          theNucleus->getStore()->getBook()->incrementBlockedDecays();
          break;
        case NoEnergyConservationFS:
        case ParticleBelowFermiFS:
        case ParticleBelowZeroFS:
          break;
        case ValidFS:
          theNucleus->getStore()->getBook()->incrementAcceptedDecays();
      }
    }
    return fs;
  }

  std::string DecayAvatar::dump() const {
    std::stringstream ss;
    ss << "(avatar " << theTime << " 'decay" << std::endl
      << "(list " << std::endl 
      << particle1->dump()
      << "))" << std::endl;
    return ss.str();
  }
}
