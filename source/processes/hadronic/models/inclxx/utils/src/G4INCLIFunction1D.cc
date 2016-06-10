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

/** \file G4INCLIFunction1D.cc
 * \brief Functor for 1-dimensional mathematical functions
 *
 * \date 16 July 2011
 * \author Davide Mancusi
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "G4INCLIFunction1D.hh"
#include "G4INCLLogger.hh"
#include "G4INCLInvFInterpolationTable.hh"

namespace G4INCL {

  const G4double IFunction1D::integrationCoefficients[] = {
      2.*95.0/288.0,
      317.0/240.0,
      23.0/30.0,
      793.0/720.0,
      157.0/160.0,
      157.0/160.0,
      793.0/720.0,
      23.0/30.0,
      317.0/240.0,
  };

  G4double IFunction1D::integrate(const G4double x0, const G4double x1, const G4double step) const {
    G4double xi = std::max(x0, xMin);
    G4double xa = std::min(x1, xMax);
    G4double sign;

    if(x1 <= x0) {
      sign = -1.0;
      std::swap(xi, xa);
    } else
      sign = 1.0;

    const G4double interval = xa - xi;

    G4int nIntervals;
    if(step<0.) {
      nIntervals = 45;
    } else {
      nIntervals = G4int(interval/step);

      // Round up nIntervals to the closest multiple of 9
      G4int remainder = nIntervals % 9;
      if (remainder != 0)
        nIntervals += 9 - remainder;

      nIntervals = std::max(nIntervals, 9);
    }

    const G4double dx = interval/nIntervals;
    G4double result = (operator()(xi) + operator()(xa)) * integrationCoefficients[0]/2;
    for(G4int j = 1; j<nIntervals; ++j) {
      const G4double x = xi + interval*G4double(j)/G4double(nIntervals);
      const unsigned index = j%9;
      result += operator()(x) * integrationCoefficients[index];
    }

    return result*dx*sign;

  }

  IFunction1D *IFunction1D::primitive() const {
    class Primitive : public IFunction1D {
      public:
        Primitive(IFunction1D const * const f) :
          IFunction1D(f->getXMinimum(), f->getXMaximum()),
          theFunction(f)
      {}

        G4double operator()(const G4double x) const {
          return theFunction->integrate(xMin,x);
        }
      private:
        IFunction1D const * const theFunction;
    } *thePrimitive = new Primitive(this);

    return thePrimitive;
  }

  InterpolationTable *IFunction1D::inverseCDFTable(IFunction1D::ManipulatorFunc fWrap, const G4int nNodes) const {
    class InverseCDF : public IFunction1D {
      public:
        InverseCDF(IFunction1D const * const f, ManipulatorFunc fw) :
          IFunction1D(f->getXMinimum(), f->getXMaximum()),
          theFunction(f),
          normalisation(1./theFunction->integrate(xMin,xMax)),
          fWrap(fw)
      {}

        G4double operator()(const G4double x) const {
          if(fWrap)
            return fWrap(std::min(1., normalisation * theFunction->integrate(xMin,x)));
          else
            return std::min(1., normalisation * theFunction->integrate(xMin,x));
        }
      private:
        IFunction1D const * const theFunction;
        const G4double normalisation;
        ManipulatorFunc fWrap;
    } *theInverseCDF = new InverseCDF(this, fWrap);

    InterpolationTable *theTable = new InvFInterpolationTable(*theInverseCDF, nNodes);
    delete theInverseCDF;
    return theTable;
  }
}

