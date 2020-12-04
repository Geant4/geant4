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
// G4TClassicalRK4
//
// Class description:
//
// Templated version of G4ClassicalRK4
//
//
// Created: Josh Xie  (supported by Google Summer of Code 2014 )
// Supervisors:  Sandro Wenzel, John Apostolakis (CERN)
// Adapted from G4G4TClassicalRK4 class
// --------------------------------------------------------------------
#include "G4ThreeVector.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4TMagErrorStepper.hh"

template <class T_Equation, unsigned int N>
class G4TClassicalRK4
  : public G4TMagErrorStepper<G4TClassicalRK4<T_Equation, N>, T_Equation, N>
{
 public:  // with description
  static constexpr G4double IntegratorCorrection = 1. / ((1 << 4) - 1);

  G4TClassicalRK4(T_Equation* EqRhs, G4int numberOfVariables = 8);

  virtual ~G4TClassicalRK4() { ; }

  void RightHandSideInl(G4double y[],
                     G4double dydx[])
  {
    fEquation_Rhs->T_Equation::RightHandSide(y, dydx);
  }

  // A stepper that does not know about errors.
  // It is used by the MagErrorStepper stepper.

  inline  // __attribute__((always_inline))
  void DumbStepper(const G4double yIn[],
                   const G4double dydx[],
                   G4double h,
                   G4double yOut[]);

 public:  // without description
  G4int IntegratorOrder() const { return 4; }

  G4TClassicalRK4(const G4TClassicalRK4&) = delete;
  G4TClassicalRK4& operator=(const G4TClassicalRK4&) = delete;
  // No copy constructor and assignment operator.

 private:
  // G4int fNumberOfVariables ; // is set default to 6 in constructor
  G4double dydxm[N < 8 ? 8 : N];
  G4double dydxt[N < 8 ? 8 : N];
  G4double yt[N < 8 ? 8 : N];
  // scratch space - not state

  T_Equation* fEquation_Rhs;
};

template <class T_Equation, unsigned int N >
G4TClassicalRK4<T_Equation,N>::
G4TClassicalRK4(T_Equation* EqRhs, G4int numberOfVariables)
   : G4TMagErrorStepper<G4TClassicalRK4<T_Equation, N>, T_Equation, N>(
      EqRhs, numberOfVariables > 8 ? numberOfVariables : 8 )
   , fEquation_Rhs(EqRhs)
{
  // unsigned int noVariables = std::max(numberOfVariables, 8);  // For Time .. 7+1
  if( dynamic_cast<G4EquationOfMotion*>(EqRhs) == nullptr )
  {
    G4Exception("G4TClassicalRK4: constructor", "GeomField0001",
                FatalException, "Equation is not an G4EquationOfMotion.");      
  }
}

template <class T_Equation, unsigned int N >
void
G4TClassicalRK4<T_Equation,N>::DumbStepper(const G4double yIn[],
                                           const G4double dydx[],
                                           G4double h,
                                           G4double yOut[]) 
// Given values for the variables y[0,..,n-1] and their derivatives
// dydx[0,...,n-1] known at x, use the classical 4th Runge-Kutta
// method to advance the solution over an interval h and return the
// incremented variables as yout[0,...,n-1], which not be a distinct
// array from y. The user supplies the routine RightHandSide(x,y,dydx),
// which returns derivatives dydx at x. The source is routine rk4 from
// NRC p. 712-713 .
{
  G4double hh = h * 0.5, h6 = h / 6.0;
  
  // Initialise time to t0, needed when it is not updated by the integration.
  //        [ Note: Only for time dependent fields (usually electric)
  //                  is it neccessary to integrate the time.]
  yt[7]   = yIn[7];
  yOut[7] = yIn[7];
  
  for(unsigned int i = 0; i < N; ++i)
  {
     yt[i] = yIn[i] + hh * dydx[i];  // 1st Step K1=h*dydx
  }
  this->RightHandSideInl(yt, dydxt);  // 2nd Step K2=h*dydxt
  
  for(unsigned int i = 0; i < N; ++i)
  {
     yt[i] = yIn[i] + hh * dydxt[i];
  }
  this->RightHandSideInl(yt, dydxm);  // 3rd Step K3=h*dydxm
  
  for(unsigned int i = 0; i < N; ++i)
  {
     yt[i] = yIn[i] + h * dydxm[i];
     dydxm[i] += dydxt[i];  // now dydxm=(K2+K3)/h
  }
  this->RightHandSideInl(yt, dydxt);  // 4th Step K4=h*dydxt
  
  for(unsigned int i = 0; i < N; ++i)  // Final RK4 output
  {
     yOut[i] = yIn[i] + h6 * (dydx[i] + dydxt[i] +
                              2.0 * dydxm[i]);  //+K1/6+K4/6+(K2+K3)/3
  }
  if(N == 12)
  {
     this->NormalisePolarizationVector(yOut);
  }
  
}  // end of DumbStepper ....................................................

// template <class T_Equation, unsigned int N >
// G4TClassicalRK4<T_Equation,N>::

