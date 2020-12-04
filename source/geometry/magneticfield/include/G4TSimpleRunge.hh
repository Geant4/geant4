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
// G4TSimpleRunge
//
// Class description:
//
// Templated version of G4SimpleRunge
//
//
// Created: Josh Xie  (supported by Google Summer of Code 2014 )
// Supervisors:  Sandro Wenzel, John Apostolakis (CERN)
// Adapted from G4G4TSimpleRunge class
// --------------------------------------------------------------------
#ifndef G4TSimpleRunge_HH
#define G4TSimpleRunge_HH

#include <cassert>
#include "G4TMagErrorStepper.hh"
#include "G4ThreeVector.hh"

template <class T_Equation, int N>
class G4TSimpleRunge
  : public G4TMagErrorStepper<G4TSimpleRunge<T_Equation, N>, T_Equation, N>
{
 public:  // with description
  static constexpr double IntegratorCorrection = 1. / ((1 << 2) - 1);

  G4TSimpleRunge(T_Equation* EqRhs, G4int numberOfVariables = 6)
    : G4TMagErrorStepper<G4TSimpleRunge<T_Equation, N>, T_Equation, N>(
        EqRhs, numberOfVariables)
    , fNumberOfVariables(numberOfVariables)
    , fEquation_Rhs(EqRhs)
      
  {
    // default GetNumberOfStateVariables() == 12
    assert(this->GetNumberOfStateVariables() <= 12);
  }

  ~G4TSimpleRunge() { ; }

  inline void RightHandSide(G4double y[],
                            G4double dydx[])
  {
    fEquation_Rhs->T_Equation::RightHandSide(y, dydx);
  }

  inline void DumbStepper(const G4double yIn[],
                          const G4double dydx[],
                          G4double h, G4double yOut[]) // override final
  {
    // Initialise time to t0, needed when it is not updated by the integration.
    yTemp[7] = yOut[7] = yIn[7];  //  Better to set it to NaN;  // TODO

    for(G4int i = 0; i < N; ++i)
    {
      yTemp[i] = yIn[i] + 0.5 * h * dydx[i];
    }

    this->RightHandSide(yTemp, dydxTemp);

    for(G4int i = 0; i < N; ++i)
    {
      yOut[i] = yIn[i] + h * (dydxTemp[i]);
    }
  }

 public:  // without description
  inline G4int IntegratorOrder() const { return 2; }

 private:
  G4int fNumberOfVariables;
  G4double dydxTemp[N > 12 ? N : 12];
  G4double yTemp[N > 12 ? N : 12];

  T_Equation* fEquation_Rhs;
  // scratch space
};

#endif /* G4TSimpleRunge_HH */
