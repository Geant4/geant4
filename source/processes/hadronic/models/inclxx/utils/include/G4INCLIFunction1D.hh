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

/** \file G4INCLIFunction1D.hh
 * \brief Functor for 1-dimensional mathematical functions
 *
 * \date 16 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLIFUNCTION1D_HH_
#define G4INCLIFUNCTION1D_HH_ 1

#include <vector>

namespace G4INCL {

  // Forward declaration
  class InterpolationTable;

  /**
   * 1D function interface
   */
  class IFunction1D {
    public:
      IFunction1D() :
        xMin(0.),
        xMax(0.)
    {};
      IFunction1D(const G4double x0, const G4double x1) :
        xMin(x0),
        xMax(x1)
    {};

      virtual ~IFunction1D() {};

      /// \brief Return the minimum allowed value of the independent variable
      virtual inline G4double getXMinimum() const { return xMin; }

      /// \brief Return the maximum allowed value of the independent variable
      virtual inline G4double getXMaximum() const { return xMax; }

      /// \brief Compute the value of the function
      virtual G4double operator()(const G4double x) const = 0;

      /** \brief Integrate the function between two values
       *
       * \param x0 lower integration bound
       * \param x1 upper integration bound
       * \param step largest integration step size; if <0, 45 steps will be used
       * \return \f$\int_{x_0}^{x_1} f(x) dx\f$
       */
      virtual G4double integrate(const G4double x0, const G4double x1, const G4double step=-1.) const;

      /// \brief Return a pointer to the (numerical) primitive to this function
      IFunction1D *primitive() const;

      /// \brief Typedef to simplify the syntax of inverseCDFTable
      typedef G4double (* const ManipulatorFunc)(const G4double);

      /** \brief Return a pointer to the inverse of the CDF of this function
       *
       * The function parameter fWrap is wrapped around the return value of
       * operator(). If fWrap=NULL (default), fWrap=identity.
       */
      InterpolationTable *inverseCDFTable(ManipulatorFunc fWrap=0, const G4int nNodes=60) const;

    protected:
      /// \brief Minimum value of the independent variable
      G4double xMin;
      /// \brief Maximum value of the independent variable
      G4double xMax;

    private:
      /// \brief Coefficients for numerical integration
      static const G4double integrationCoefficients[];
  };

}

#endif // G4INCLIFUNCTION1D_HH_
