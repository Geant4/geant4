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

/** \file G4INCLCrossSectionsTruncatedMultiPions.hh
 * \brief Truncated multipion cross sections
 *
 * \date 10th December 2014
 * \author Davide Mancusi
 */

#ifndef G4INCLCROSSSECTIONSTRUNCATEDMULTIPIONS_HH
#define G4INCLCROSSSECTIONSTRUNCATEDMULTIPIONS_HH

#include "G4INCLCrossSectionsMultiPions.hh"
#include <limits>

namespace G4INCL {
  /// \brief Truncated multipion cross sections

  class CrossSectionsTruncatedMultiPions : public CrossSectionsMultiPions {
    public:
      CrossSectionsTruncatedMultiPions(const G4int nPi=std::numeric_limits<G4int>::max());

      /** \brief Elastic particle-particle cross section
       *
       * The elastic pi-N cross section is calculated by difference:
       *   elastic = total - inelastic(pion production) - delta,
       * but we need to use the unmodified piNToDelta cross section, or else we
       * will violate the cross-section sum.
       *
       * This modification is only necessary for piN collisions.
       */
      virtual G4double elastic(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for Delta production - piN Channel
      virtual G4double piNToDelta(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for X pion production - piN Channel
      virtual G4double piNToxPiN(const G4int xpi, Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for X pion production - NN Channel
      virtual G4double NNToxPiNN(const G4int xpi, Particle const * const p1, Particle const * const p2);

    protected:

      /// \brief Maximum number of pions produced in TruncatedMultiPion collisions
      const G4int nMaxPi;

  };
}

#endif
