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
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4INCLIFunction_hh
#define G4INCLIFunction_hh 1

namespace G4INCL {

  /**
   * 1D function G4interface
   */
  class IFunction1D {
  public:
    IFunction1D() {};
    IFunction1D(G4double, G4double, G4double) {};
    virtual ~IFunction1D() {};

    /**
     * Compute the value of the function
     */
    virtual G4double getValue(G4double r) = 0;

    virtual G4double getRadiusParameter() = 0;
    virtual G4double getMaximumRadius() = 0;
    virtual G4double getDiffusenessParameter() = 0;
    virtual void setRadiusParameter(G4double r) = 0;
    virtual void setMaximumRadius(G4double r) = 0;
    virtual void setDiffusenessParameter(G4double a) = 0;
  };

  class WoodsSaxon : public IFunction1D {
  public:
    WoodsSaxon(G4double radiusParameter, G4double maximumRadius, G4double diffusenessParameter)
      :theRadiusParameter(radiusParameter), theMaximumRadius(maximumRadius),
       theDiffusenessParameter(diffusenessParameter)
    {};

    virtual ~WoodsSaxon() {};

    /**
     * r^2 / (1.0 - exp((r-r0)/adif))
     */
    inline G4double getValue(G4double r) {
      return theRadiusParameter*theRadiusParameter
        / (1.0 + std::exp((r - theRadiusParameter)/theDiffusenessParameter));
    };

    inline G4double getRadiusParameter() { return theRadiusParameter; };
    inline G4double getMaximumRadius() { return theMaximumRadius; };
    inline G4double getDiffusenessParameter() { return theDiffusenessParameter; };

    inline void setRadiusParameter(G4double r) { theRadiusParameter = r; };
    inline void setMaximumRadius(G4double r) { theMaximumRadius = r; };
    inline void setDiffusenessParameter(G4double a) { theDiffusenessParameter = a; };

  private:
    G4double theRadiusParameter, theMaximumRadius, theDiffusenessParameter;
  };

  class DerivWoodsSaxon : public IFunction1D {
  public:
    DerivWoodsSaxon(G4double radiusParameter, G4double maximumRadius, G4double diffusenessParameter)
      :theRadiusParameter(radiusParameter), theMaximumRadius(maximumRadius),
       theDiffusenessParameter(diffusenessParameter)
    {};

    virtual ~DerivWoodsSaxon() {};

    inline G4double getValue(G4double r) {
      G4double derivwsax = std::pow(r,3)
	*std::exp((r - theRadiusParameter)/theDiffusenessParameter)
	/std::pow((1.0 + std::exp((r - theRadiusParameter)/theDiffusenessParameter)),2);
      return derivwsax/theDiffusenessParameter;
    }

    inline G4double getRadiusParameter() { return theRadiusParameter; };
    inline G4double getMaximumRadius() { return theMaximumRadius; };
    inline G4double getDiffusenessParameter() { return theDiffusenessParameter; };

    inline void setRadiusParameter(G4double r) { theRadiusParameter = r; };
    inline void setMaximumRadius(G4double r) { theMaximumRadius = r; };
    inline void setDiffusenessParameter(G4double a) { theDiffusenessParameter = a; };

  private:
    G4double theRadiusParameter, theMaximumRadius, theDiffusenessParameter;
  };

  class DerivModifiedHarmonicOscillator : public IFunction1D {
  public:
    DerivModifiedHarmonicOscillator(G4double radiusParameter, G4double maximumRadius, G4double diffusenessParameter)
      :theRadiusParameter(radiusParameter), theMaximumRadius(maximumRadius),
       theDiffusenessParameter(diffusenessParameter)
    {};

    virtual ~DerivModifiedHarmonicOscillator() {};

    inline G4double getValue(G4double r) {
      const G4double arg = std::pow((r/theDiffusenessParameter),2);
      return -2.0* r*r *arg * (theRadiusParameter - 1.0 - theRadiusParameter*arg)*std::exp(-arg);
    }

    inline G4double getRadiusParameter() { return theRadiusParameter; };
    inline G4double getMaximumRadius() { return theMaximumRadius; };
    inline G4double getDiffusenessParameter() { return theDiffusenessParameter; };

    inline void setRadiusParameter(G4double r) { theRadiusParameter = r; };
    inline void setMaximumRadius(G4double r) { theMaximumRadius = r; };
    inline void setDiffusenessParameter(G4double a) { theDiffusenessParameter = a; };

  private:
    G4double theRadiusParameter, theMaximumRadius, theDiffusenessParameter;
  };

  class DerivGaussian : public IFunction1D {
  public:
    DerivGaussian(G4double radiusParameter, G4double maximumRadius, G4double diffusenessParameter)
      :theRadiusParameter(radiusParameter), theMaximumRadius(maximumRadius),
       theDiffusenessParameter(diffusenessParameter)
    {};

    virtual ~DerivGaussian() {};

    inline G4double getValue(G4double r) {
      const G4double arg = std::pow((r/theDiffusenessParameter),2);
      return r*r*arg*std::exp(-arg/2.0);
    }

    inline G4double getRadiusParameter() { return theRadiusParameter; };
    inline G4double getMaximumRadius() { return theMaximumRadius; };
    inline G4double getDiffusenessParameter() { return theDiffusenessParameter; };

    inline void setRadiusParameter(G4double r) { theRadiusParameter = r; };
    inline void setMaximumRadius(G4double r) { theMaximumRadius = r; };
    inline void setDiffusenessParameter(G4double a) { theDiffusenessParameter = a; };

  private:
    G4double theRadiusParameter, theMaximumRadius, theDiffusenessParameter;
  };
}

#endif
