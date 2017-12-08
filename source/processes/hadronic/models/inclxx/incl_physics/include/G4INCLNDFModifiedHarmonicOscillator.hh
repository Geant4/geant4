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

/** \file G4INCLNDFModifiedHarmonicOscillator.hh
 * \brief Class for modified harmonic oscillator density
 *
 * \date 16 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLNDFMODIFIEDHARMONICOSCILLATOR_HH_
#define G4INCLNDFMODIFIEDHARMONICOSCILLATOR_HH_

#include "G4INCLIFunction1D.hh"
#include <cmath>
#include <algorithm>

namespace G4INCL {

  namespace NuclearDensityFunctions {

    class ModifiedHarmonicOscillatorRP : public IFunction1D {
      public:
        ModifiedHarmonicOscillatorRP(G4double radiusParameter, G4double maximumRadius, G4double diffusenessParameter) :
          IFunction1D(0., maximumRadius),
          theRadiusParameter(radiusParameter),
          theDiffusenessParameter(diffusenessParameter)
      {}

        inline G4double operator()(const G4double r) const {
          const G4double arg = std::pow((r/theDiffusenessParameter),2);
          return std::max(0., -2.0* r*r *arg * (theRadiusParameter - 1.0 - theRadiusParameter*arg)*std::exp(-arg));
        }

        inline G4double getRadiusParameter() { return theRadiusParameter; };
        inline G4double getDiffusenessParameter() { return theDiffusenessParameter; };

        inline void setRadiusParameter(G4double r) { theRadiusParameter = r; };
        inline void setDiffusenessParameter(G4double a) { theDiffusenessParameter = a; };

      protected:
        G4double theRadiusParameter, theDiffusenessParameter;
    };

    class ModifiedHarmonicOscillator : public IFunction1D {
      public:
        ModifiedHarmonicOscillator(G4double radiusParameter, G4double maximumRadius, G4double diffusenessParameter) :
          IFunction1D(0., maximumRadius),
          theRadiusParameter(radiusParameter),
          theDiffusenessParameter(diffusenessParameter),
          normalisation(2./((theDiffusenessParameter+theRadiusParameter)*std::pow(theDiffusenessParameter,2.)))
      {}

        inline G4double operator()(const G4double r) const {
          const G4double arg = std::pow((r/theDiffusenessParameter),2);
          return r*r*( 1. + theRadiusParameter * arg) * std::exp(-arg);
        }

        inline G4double getRadiusParameter() { return theRadiusParameter; };
        inline G4double getDiffusenessParameter() { return theDiffusenessParameter; };

        inline void setRadiusParameter(G4double r) { theRadiusParameter = r; };
        inline void setDiffusenessParameter(G4double a) { theDiffusenessParameter = a; };

      protected:
        G4double theRadiusParameter, theDiffusenessParameter;
        const G4double normalisation;
    };

  }

}

#endif // G4INCLNDFMODIFIEDHARMONICOSCILLATOR_HH_

