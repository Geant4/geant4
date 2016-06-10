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

/** \file G4INCLDeuteronDensity.cc
 * \brief Deuteron density in r and p according to the Paris potential.
 *
 * \date 6 March 2012
 * \author Davide Mancusi
 */

#include "G4INCLDeuteronDensity.hh"
#include "G4INCLGlobals.hh"
// #include <cassert>
#include <algorithm>

namespace G4INCL {

  namespace DeuteronDensity {

    namespace {

      const G4int coeffTableSize = 13;

      /// \brief Coefficients for the deuteron wave function
      const G4double coeff1[coeffTableSize] = {
        0.88688076e+0,
        -0.34717093e+0,
        -0.30502380e+1,
        0.56207766e+2,
        -0.74957334e+3,
        0.53365279e+4,
        -0.22706863e+5,
        0.60434469e+5,
        -0.10292058e+6,
        0.11223357e+6,
        -0.75925226e+5,
        0.29059715e+5,
        -0.48157368e+4
      };

      /// \brief Coefficients for the deuteron wave function
      const G4double coeff2[coeffTableSize] = {
        0.23135193e-1,
        -0.85604572e+0,
        0.56068193e+1,
        -0.69462922e+2,
        0.41631118e+3,
        -0.12546621e+4,
        0.12387830e+4,
        0.33739172e+4,
        -0.13041151e+5,
        0.19512524e+5,
        -0.15634324e+5,
        0.66231089e+4,
        -0.11698185e+4
      };

      /// \brief Normalisation coefficient for the r-space deuteron wave function
      const G4double normalisationR = std::sqrt(32. * Math::pi) * 0.28212;

      /// \brief Normalisation coefficient for the p-space deuteron wave function
      const G4double normalisationP = normalisationR / (std::sqrt(4. * Math::pi) * std::pow(PhysicalConstants::hc,1.5));

      /// \brief Mysterious coefficient that appears in the wavefunctions
      const G4double al = 0.23162461;

    }

    G4double densityR(const G4double r) {
      const G4double sWave = wavefunctionR(0, r);
      const G4double dWave = wavefunctionR(2, r);
      return r*r*(sWave*sWave + dWave*dWave);
    }

    G4double derivDensityR(const G4double r) {
      const G4double sWave = wavefunctionR(0, r);
      const G4double dWave = wavefunctionR(2, r);
      const G4double sWaveDeriv = derivWavefunctionR(0, r);
      const G4double dWaveDeriv = derivWavefunctionR(2, r);
      return (sWave*sWaveDeriv + dWave*dWaveDeriv) / Math::twoPi;
    }

    G4double densityP(const G4double p) {
      const G4double sWave = wavefunctionP(0, p);
      const G4double dWave = wavefunctionP(2, p);
      return p*p*(sWave*sWave + dWave*dWave);
    }

    G4double wavefunctionR(const G4int l, const G4double theR) {
// assert(l==0 || l==2); // only s- and d-waves in a deuteron
      const G4double r = 2. * std::max(theR, 1.e-4);

      G4double result = 0.;
      G4double fmr;

      for(G4int i=0; i<coeffTableSize; ++i) {
        fmr = r * (al+i);
        if(l==0) { // s-wave
          result += coeff1[i] * std::exp(-fmr);
        } else { // d-wave
          result += coeff2[i] * std::exp(-fmr) * (1.+3./fmr+3./(fmr*fmr));
        }
      }

      result *= normalisationR/r;
      return result;
    }

    G4double derivWavefunctionR(const G4int l, const G4double theR) {
// assert(l==0 || l==2); // only s- and d-waves in a deuteron
      const G4double r = 2. * std::max(theR, 1.e-4);

      G4double result = 0.;
      G4double fmr;

      for(G4int i=0; i<coeffTableSize; ++i) {
        fmr = r * (al+i);
        if(l==0) { // s-wave
          result += coeff1[i] * std::exp(-fmr) * (fmr + 1.);
        } else { // d-wave
          result += coeff2[i] * std::exp(-fmr) * (fmr + 4. + 9./fmr + 9./(fmr*fmr));
        }
      }

      result *= -normalisationR/(r*r);
      return result;
    }

    G4double wavefunctionP(const G4int l, const G4double theQ) {
// assert(l==0 || l==2); // only s- and d-waves in a deuteron
      const G4double q = theQ / PhysicalConstants::hc;
      const G4double q2 = q*q;
      G4double result=0.;
      G4double fmq, alPlusI;
      for(G4int i=0; i<coeffTableSize; ++i) {
        alPlusI = al+i;
        fmq = q2 + alPlusI*alPlusI;
        if(l==0) { // s-wave
          result += coeff1[i] / fmq;
        } else { // d-wave
          result += coeff2[i] / fmq;
        }
      }

      result *= normalisationP;
      return result;
    }

  }

}
