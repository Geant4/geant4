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

/** \file G4INCLRootFinder.hh
 * \brief Static root-finder algorithm.
 *
 * Provides a (nearly) stateless root-finder algorithm.
 *
 * \date 2nd March 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLROOTFINDER_HH_
#define G4INCLROOTFINDER_HH_

#include <utility>
#include "G4INCLIFunction1D.hh"

namespace G4INCL {

  class RootFunctor : public IFunction1D {
  public:
    virtual void cleanUp(const G4bool success) const = 0;
    virtual ~RootFunctor() {};
  protected:
    RootFunctor(const G4double x0, const G4double x1) :
      IFunction1D(x0, x1)
    {};
  };

  namespace RootFinder {

    class Solution {
      public:
        Solution() :
          success(false),
          x(0.),
          y(0.)
      {}
        Solution( const G4double x0, const G4double y0) :
          success(true),
          x(x0),
          y(y0)
      {}
        ~Solution() {}

        G4bool success;
        G4double x;
        G4double y;
    };

    /** \brief Numerically solve a one-dimensional equation.
     *
     * Numerically solves the equation f(x)==0. This implementation uses the
     * false-position method.
     *
     * If a root is found, it can be retrieved using the getSolution() method,
     *
     * \param f pointer to a RootFunctor
     * \param x0 initial value of the function argument
     * \return a Solution object describing the root, if it was found
     */
    Solution solve(RootFunctor const * const f, const G4double x0);

  }
}
#endif /* G4INCLROOTFINDER_HH_ */
