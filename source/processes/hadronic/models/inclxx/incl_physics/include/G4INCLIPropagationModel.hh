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
 * IPropagationModel.hh
 *
 *  \date 4 juin 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLIPropagationModel_hh
#define G4INCLIPropagationModel_hh

#include "G4INCLIAvatar.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLFinalState.hh"

namespace G4INCL {

    /**
     * Propagation model takes care of transporting the particles until something
     * interesting (i.e. an avatar) happens. This avatar is then returned back to the INCL
     * kernel for further processing.
     *
     * The propagation model idea abstracts the details of propagation. This allows us to
     * conveniently support multiple propagation models and to compare their results. Some
     * possible future propagation models are: straight line trajectories by using constant
     * time step and curved trajectories.
     */
    class IPropagationModel {
    public:
      IPropagationModel() {}
      virtual ~IPropagationModel() {}

      /**
       * Set the nucleus for the propagation model.
       *
       * @param nucleus Pointer to the nucleus
       */
      virtual void setNucleus(G4INCL::Nucleus *nucleus) = 0;

      /**
       * Get a pointer to the nucleus.
       *
       * @return G4INCL::Nuleus*
       */
      virtual G4INCL::Nucleus* getNucleus() = 0;

      virtual G4double shoot(ParticleSpecies const &projectileSpecies, const G4double kineticEnergy, const G4double impactParameter, const G4double phi) = 0;
    protected:
      virtual G4double shootParticle(ParticleType const t, const G4double kineticEnergy, const G4double impactParameter, const G4double phi) = 0;
      virtual G4double shootComposite(ParticleSpecies const &s, const G4double kineticEnergy, const G4double impactParameter, const G4double phi) = 0;

    public:

      /**
       * Returns the current global time of the system.
       */
      virtual G4double getCurrentTime() = 0;

      /**
       * Set new stopping time to the propagation.
       */
      virtual void setStoppingTime(G4double) = 0;

      /**
       * Get the current stopping time.
       */
      virtual G4double getStoppingTime() = 0;

      /**
       * Propagate the particles and get the next avatar.
       *
       * @return G4INCL::IAvatar the next avatar
       */
      virtual G4INCL::IAvatar* propagate(FinalState const * const fs) = 0;
    };

}

#endif
