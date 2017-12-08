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

#include "G4INCLDecayAvatar.hh"

#include "G4INCLDeltaDecayChannel.hh"
#include "G4INCLPionResonanceDecayChannel.hh"
#include "G4INCLSigmaZeroDecayChannel.hh"
#include "G4INCLNeutralKaonDecayChannel.hh"
#include "G4INCLStrangeAbsorbtionChannel.hh"
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
  
  DecayAvatar::DecayAvatar(G4INCL::Particle *aParticle, G4INCL::Particle *bParticle, G4double time, G4INCL::Nucleus *n, G4bool force)
    : InteractionAvatar(time, n, aParticle, bParticle), forced(force),
      incidentDirection(aParticle->getMomentum())
  {
    setType(DecayAvatarType);
  }

  DecayAvatar::~DecayAvatar() {

  }

  G4INCL::IChannel* DecayAvatar::getChannel() {
    if(!particle2){
       if(particle1->isDelta()) {
         INCL_DEBUG("DeltaDecayChannel chosen." << '\n');
         return new DeltaDecayChannel(particle1, incidentDirection);
       }
       else if(particle1->isEta() || particle1->isOmega()) {
         INCL_DEBUG("PionResonanceDecayChannel chosen." << '\n');
         return new PionResonanceDecayChannel(particle1, incidentDirection);
       }
       else if(particle1->getType() == SigmaZero) {
         INCL_DEBUG("SigmaZeroDecayChannel chosen." << '\n');
         return new SigmaZeroDecayChannel(particle1, incidentDirection);
       }
       else if(particle1->getType() == KZero || particle1->getType() == KZeroBar) {
         INCL_DEBUG("NeutralKaonDecayChannel chosen." << '\n');
         return new NeutralKaonDecayChannel(particle1);
       }
    }
    else if(((particle1->isAntiKaon() || particle1->isSigma()) && particle2->isNucleon()) || ((particle2->isAntiKaon() || particle2->isSigma()) && particle1->isNucleon())){
      INCL_DEBUG("StrangeAbsorbtion." << '\n');
      return new StrangeAbsorbtionChannel(particle1, particle2);
    }
    return NULL;
  }

  void DecayAvatar::preInteraction() {
    InteractionAvatar::preInteraction();
  }

  void DecayAvatar::postInteraction(FinalState *fs) {
    // Make sure we have at least two particles in the final state
    // Removed because of neutral kaon decay
    
// assert((fs->getModifiedParticles().size() + fs->getCreatedParticles().size() - fs->getDestroyedParticles().size() >= 2) || ((*fs->getModifiedParticles().begin())->getType() == KShort || (*fs->getModifiedParticles().begin())->getType() == KLong ));
    //assert((fs->getModifiedParticles().size() + fs->getCreatedParticles().size() - fs->getDestroyedParticles().size() >= 1));
    if(!forced) { // Normal decay
      // Call the postInteraction method of the parent class
      // (provides Pauli blocking and enforces energy conservation)
      InteractionAvatar::postInteraction(fs);

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
        fs->addModifiedParticle(particle1);
    } else { // Forced decay
      modified = fs->getModifiedParticles();
      created = fs->getCreatedParticles();
      Destroyed = fs->getDestroyedParticles();
      modifiedAndCreated = modified;
      modifiedAndCreated.insert(modifiedAndCreated.end(), created.begin(), created.end());
      ModifiedAndDestroyed = modified;
      ModifiedAndDestroyed.insert(ModifiedAndDestroyed.end(), Destroyed.begin(), Destroyed.end());
      
      std::vector<G4int> newBiasCollisionVector;
      newBiasCollisionVector = ModifiedAndDestroyed.getParticleListBiasVector();
      for(ParticleIter i=modifiedAndCreated.begin(), e=modifiedAndCreated.end(); i!=e; ++i ) {
	    (*i)->setBiasCollisionVector(newBiasCollisionVector);
      }
      // Try to enforce energy conservation
      fs->setTotalEnergyBeforeInteraction(oldTotalEnergy);
      const G4bool success = enforceEnergyConservation(fs);
      if(!success) {
        INCL_DEBUG("Enforcing energy conservation: failed!" << '\n');

        if(theNucleus) {
          // Restore the state of the initial particles
          restoreParticles();

          // Delete newly created particles
          for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i )
            delete *i;

          fs->reset();
          fs->makeNoEnergyConservation();
          fs->setTotalEnergyBeforeInteraction(0.0);

          return; // Interaction is blocked. Return an empty final state.
        } else {
          // If there is no nucleus we have to continue anyway, even if energy
          // conservation failed. We cannot afford producing unphysical
          // remnants.
          INCL_DEBUG("No nucleus, continuing anyway." << '\n');
        }
      } else {
        INCL_DEBUG("Enforcing energy conservation: success!" << '\n');
      }

      if(theNucleus) {
        // Test CDPP blocking
        G4bool isCDPPBlocked = Pauli::isCDPPBlocked(created, theNucleus);

        if(isCDPPBlocked) {
          INCL_DEBUG("CDPP: Blocked!" << '\n');

          // Restore the state of both particles
          restoreParticles();

          // Delete newly created particles
          for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i )
            delete *i;

          fs->reset();
          fs->makePauliBlocked();
          fs->setTotalEnergyBeforeInteraction(0.0);

          return; // Interaction is blocked. Return an empty final state.
        }
        INCL_DEBUG("CDPP: Allowed!" << '\n');

      }
    }
    // If there is a nucleus, increment the counters
    if(theNucleus) {
      switch(fs->getValidity()) {
        case PauliBlockedFS:
          theNucleus->getStore()->getBook().incrementBlockedDecays();
          break;
        case NoEnergyConservationFS:
        case ParticleBelowFermiFS:
        case ParticleBelowZeroFS:
          break;
        case ValidFS:
          theNucleus->getStore()->getBook().incrementAcceptedDecays();
      }
    }
    
    return;
  }

  std::string DecayAvatar::dump() const {
    std::stringstream ss;
    ss << "(avatar " << theTime << " 'decay" << '\n'
      << "(list " << '\n'
      << particle1->dump()
      << "))" << '\n';
    return ss.str();
  }
}
