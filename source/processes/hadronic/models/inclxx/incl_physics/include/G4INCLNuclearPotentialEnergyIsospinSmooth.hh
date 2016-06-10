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

/** \file G4INCLNuclearPotentialEnergyIsospinSmooth.hh
 * \brief Isospin- and energy-dependent nuclear potential.
 *
 * Provides an isospin- and energy-dependent nuclear potential. The cusp at 200
 * MeV is replaced by a smooth exponential.
 *
 * \date 16 February 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLNUCLEARPOTENTIALENERGYISOSPINSMOOTH_HH
#define G4INCLNUCLEARPOTENTIALENERGYISOSPINSMOOTH_HH 1

#include "G4INCLNuclearPotentialIsospin.hh"

namespace G4INCL {

  namespace NuclearPotential {

    class NuclearPotentialEnergyIsospinSmooth : public NuclearPotentialIsospin {

      public:
        NuclearPotentialEnergyIsospinSmooth(const G4int A, const G4int Z, const G4bool pionPotential);
        virtual ~NuclearPotentialEnergyIsospinSmooth();

        virtual G4double computePotentialEnergy(const Particle * const p) const;

      private:
        /// Slope of the V(T) curve
        static const G4double alpha;

        /// Distance from the cusp where the exponential kicks in
        static const G4double deltaE;

    };

  }
}

#endif /* G4INCLNUCLEARPOTENTIALENERGYISOSPINSMOOTH_HH */
