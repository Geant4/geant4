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
// G4TExplicitEuler
//
// Class description:
//
// Templated version of G4ExplicitEuler
//
//
// Created: Josh Xie, June 2014  (supported by Google Summer of Code 2014 )
// 
// Supervisors:  Sandro Wenzel, John Apostolakis (CERN)
// Adapted from G4G4TExplicitEuler class
// -------------------------------------------------------------------
//
// Information from G4Explicit Euler:
// ----------------------------------
// Explicit Euler stepper for magnetic field: x_1 = x_0 + h * dx_0.
// The simplistic approach to solving linear differential equations.
// Take the current derivative and add it to the current position.
// 
// Created: W.Wander <wwc@mit.edu>, 12.09.1997
// -------------------------------------------------------------------

// --------------------------------------------------------------------
#ifndef G4TExplicitEuler_HH
#define G4TExplicitEuler_HH

#include "G4TMagErrorStepper.hh"
#include "G4ThreeVector.hh"

template <class T_Equation, int N>
class G4TExplicitEuler
  : public G4TMagErrorStepper<G4TExplicitEuler<T_Equation, N>, T_Equation, N>
{
 public:  // with description
  static constexpr double IntegratorCorrection = 1.;

  G4TExplicitEuler(T_Equation* EqRhs, G4int numberOfVariables = N)
    : G4TMagErrorStepper<G4TExplicitEuler<T_Equation, N>, T_Equation, N>(
        EqRhs, numberOfVariables)
    , fEquation_Rhs(EqRhs)
  {
    if( numberOfVariables != N ){
       G4ExceptionDescription msg;
       msg << "Equation has an incompatible number of variables." ;
       msg << "   template N = " << N << " equation-Nvar= "
           << numberOfVariables;
       G4Exception("G4TExplicitEuler: constructor", "GeomField0003",
                   FatalErrorInArgument, msg );
    }
  }
   
  ~G4TExplicitEuler() { ; }

  inline void DumbStepper(const G4double yIn[],
                          const G4double dydx[],
                          G4double h, G4double yOut[]) // override final
  {
    // Initialise time to t0, needed when it is not updated by the integration.
    // yOut[7] = yIn[7];   //  Better to set it to NaN;  // TODO

    for(G4int i = 0; i < N; ++i)
    {
      yOut[i] = yIn[i] + h * dydx[i];  // 1st and only Step
    }

    return;
  }

 public:  // without description
  G4int IntegratorOrder() const { return 1; }

 private:
  T_Equation* fEquation_Rhs;
};

#endif /* G4TExplicitEuler_HH */
