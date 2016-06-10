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

/** \file G4INCLDeuteronDensity.hh
 * \brief Deuteron density in r and p according to the Paris potential.
 *
 * \date 6 March 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLDEUTERONDENSITY_HH_
#define G4INCLDEUTERONDENSITY_HH_ 1

#include "G4INCLIFunction1D.hh"

namespace G4INCL {

  /// \brief Static class for the deuteron density
  namespace DeuteronDensity {

    /** \brief PDF for a nucleon in r space
     *
     * The distribution is normalised to 1.
     *
     * \param r distance from the deuteron centre [fm]
     * \return 4 * pi * r^2 * |psi|^2
     */
    G4double densityR(const G4double r);

    /** \brief First derivative of the r-space density function
     *
     * \param r distance from the deuteron centre [fm]
     * \return d|psi|^2/dr
     */
    G4double derivDensityR(const G4double r);

    /** \brief PDF for a nucleon in p space
     *
     * The distribution is normalised to 1.
     *
     * \param p nucleon momentum [MeV/c]
     * \return 4 * pi * p^2 * d|psi|^2/dp
     */
    G4double densityP(const G4double p);

    G4double wavefunctionR(const G4int l, const G4double r);
    G4double wavefunctionP(const G4int l, const G4double p);
    G4double derivWavefunctionR(const G4int l, const G4double r);

  }

}

#endif // G4INCLDEUTERONDENSITY_HH_
