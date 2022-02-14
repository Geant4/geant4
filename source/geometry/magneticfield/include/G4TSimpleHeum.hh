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
// G4TSimpleHeum
//
// Class description:
//
// Templated version of G4SimpleHeum
//
// Created: Josh Xie  (supported by Google Summer of Code 2014 )
// Supervisors:  Sandro Wenzel, John Apostolakis (CERN)
// --------------------------------------------------------------------
// Adapted from G4G4TSimpleHeum class
// Original desription:
//
// Simple Heum stepper for magnetic field:
//        x_1 = x_0  +
//              h * 1/4 * dx(t0,x0)  +
//                  3/4 * dx(t0+2/3*h, x0+2/3*h*(dx(t0+h/3,x0+h/3*dx(t0,x0)))) 
//
// Third order solver.
// Created: W.Wander <wwc@mit.edu>, 12/09/1997
// --------------------------------------------------------------------

#ifndef TSIMPLEHEUM_HH
#define TSIMPLEHEUM_HH

#include <cassert>
#include "G4TMagErrorStepper.hh"
#include "G4ThreeVector.hh"

template <class T_Equation, unsigned int N>
class G4TSimpleHeum
  : public G4TMagErrorStepper<G4TSimpleHeum<T_Equation, N>, T_Equation, N>
{
 public:  // with description
  constexpr static unsigned int gIntegratorOrder = 3;
  static constexpr double IntegratorCorrection= 1.0 /
                                                ((1<<gIntegratorOrder) - 1);

  G4TSimpleHeum(T_Equation* EqRhs, unsigned int numberOfVariables = 6);

  ~G4TSimpleHeum() { ; }
  // Constructor and destructor.

  inline void RightHandSide(G4double y[],
                            G4double dydx[])
  {
    fEquation_Rhs->T_Equation::RightHandSide(y, dydx);
  }

  inline void DumbStepper(const G4double yIn[],
                          const G4double dydx[],
                          G4double h, G4double yOut[]); // override final

 public:  // without description
  G4int IntegratorOrder() const { return gIntegratorOrder; }

 private:
  G4int fNumberOfVariables;

  G4double dydxTemp[N];
  G4double dydxTemp2[N];
  G4double yTemp[N];
  G4double yTemp2[N];
  // scratch space
   
  T_Equation* fEquation_Rhs;
};

template <class T_Equation, unsigned int N >
G4TSimpleHeum<T_Equation,N>::G4TSimpleHeum(T_Equation* EqRhs,
                unsigned int numberOfVariables )
    : G4TMagErrorStepper<G4TSimpleHeum<T_Equation, N>, T_Equation, N>(
        EqRhs, numberOfVariables)
    , fNumberOfVariables(numberOfVariables)
    , fEquation_Rhs(EqRhs)
{
  assert(fNumberOfVariables == N);
  if( dynamic_cast<G4EquationOfMotion*>(EqRhs) == nullptr )
  {
    G4Exception("G4TSimpleHeum: constructor", "GeomField0001",
                FatalException, "Equation is not an G4EquationOfMotion.");      
  }    
}

template <class T_Equation, unsigned int N >
inline void
G4TSimpleHeum<T_Equation,N>::DumbStepper(const G4double yIn[],
                                         const G4double dydx[],
                                         G4double h, G4double yOut[])
{
  for(unsigned int i = 0; i < N; ++i)
  {
     yTemp[i] = yIn[i] + (1.0 / 3.0) * h * dydx[i];
  }
  
  this->RightHandSide(yTemp, dydxTemp);
  
  for(unsigned int i = 0; i < N; ++i)
  {
     yTemp2[i] = yIn[i] + (2.0 / 3.0) * h * dydxTemp[i];
  }
  
  this->RightHandSide(yTemp2, dydxTemp2);
  
  for(unsigned int i = 0; i < N; ++i)
  {
     yOut[i] = yIn[i] + h * (0.25 * dydx[i] + 0.75 * dydxTemp2[i]);
  }
  
  if(fNumberOfVariables == 12)
  {
     this->NormalisePolarizationVector(yOut);
  }
}

#endif /* TSIMPLEHEUM_HH */
