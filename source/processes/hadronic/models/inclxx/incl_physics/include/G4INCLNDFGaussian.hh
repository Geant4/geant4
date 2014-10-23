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

/** \file G4INCLNDFGaussian.hh
 * \brief Class for Gaussian density
 *
 * \date 16 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLNDFGAUSSIAN_HH_
#define G4INCLNDFGAUSSIAN_HH_

#include "G4INCLIFunction1D.hh"
#include <cmath>

namespace G4INCL {

  namespace NuclearDensityFunctions {

    class GaussianRP : public IFunction1D {
      public:
        GaussianRP(G4double maximumRadius, G4double standardDeviation) :
          IFunction1D(0., maximumRadius),
          theStandardDeviation(standardDeviation)
      {}

        inline G4double operator()(const G4double r) const {
          const G4double arg = std::pow((r/theStandardDeviation),2);
          return r*r*arg*std::exp(-arg/2.0);
        }

      protected:
        G4double theStandardDeviation;
    };

    class Gaussian : public IFunction1D {
      public:
        Gaussian(G4double maximumRadius, G4double standardDeviation) :
          IFunction1D(0., maximumRadius),
          theStandardDeviation(standardDeviation),
          normalisation(std::sqrt(2./Math::pi)/theStandardDeviation)
      {}

        inline G4double operator()(const G4double r) const {
          const G4double arg = std::pow((r/theStandardDeviation),2);
          return normalisation * arg * std::exp(-arg/2.0);
        }

      protected:
        G4double theStandardDeviation;
        const G4double normalisation;
    };

  }

}

#endif // G4INCLNDFGAUSSIAN_HH_

