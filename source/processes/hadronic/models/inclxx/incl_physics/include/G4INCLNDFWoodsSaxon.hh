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

/** \file G4INCLNDFWoodsSaxon.hh
 * \brief Class for Woods-Saxon density
 *
 * \date 16 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLNDFWOODSSAXON_HH_
#define G4INCLNDFWOODSSAXON_HH_

#include "G4INCLIFunction1D.hh"
#include <cmath>

namespace G4INCL {

  namespace NuclearDensityFunctions {

    class WoodsSaxonRP : public IFunction1D {
      public:
        WoodsSaxonRP(G4double radiusParameter, G4double maximumRadius, G4double diffusenessParameter) :
          IFunction1D(0., maximumRadius),
          theRadiusParameter(radiusParameter),
          theDiffusenessParameter(diffusenessParameter)
      {}

        inline G4double operator()(const G4double r) const {
          G4double wsax = std::pow(r,3)
            *std::exp((r - theRadiusParameter)/theDiffusenessParameter)
            /std::pow((1.0 + std::exp((r - theRadiusParameter)/theDiffusenessParameter)),2);
          return wsax/theDiffusenessParameter;
        }

        inline G4double getRadiusParameter() { return theRadiusParameter; };
        inline G4double getDiffusenessParameter() { return theDiffusenessParameter; };

        inline void setRadiusParameter(G4double r) { theRadiusParameter = r; };
        inline void setDiffusenessParameter(G4double a) { theDiffusenessParameter = a; };

      protected:
        G4double theRadiusParameter, theDiffusenessParameter;
    };

    class WoodsSaxon : public IFunction1D {
      public:
        WoodsSaxon(G4double radiusParameter, G4double maximumRadius, G4double diffusenessParameter) :
          IFunction1D(0., maximumRadius),
          theRadiusParameter(radiusParameter),
          theDiffusenessParameter(diffusenessParameter)
      {}

        inline G4double operator()(const G4double r) const {
          return r * r / (1.0 + std::exp((r - theRadiusParameter)/theDiffusenessParameter));
        }

        inline G4double getRadiusParameter() { return theRadiusParameter; };
        inline G4double getDiffusenessParameter() { return theDiffusenessParameter; };

        inline void setRadiusParameter(G4double r) { theRadiusParameter = r; };
        inline void setDiffusenessParameter(G4double a) { theDiffusenessParameter = a; };

      protected:
        G4double theRadiusParameter, theDiffusenessParameter;
    };

  }

}

#endif // G4INCLNDFWOODSSAXON_HH_
