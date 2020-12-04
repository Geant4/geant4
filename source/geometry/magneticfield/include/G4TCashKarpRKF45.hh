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
// G4TCashKarpRKF45 
//
// Class description:
//
// Templated version of Cash-Karp 4th/5th order embedded stepper
//
// Knowing the type (class) of the equation of motion enables a non-
// virtual call of its methods. 
// As an embedded 5th order method, it requires fewer field evaluations
// (1 initial + 5 others per step = 6 per step) than ClassicalRK4 and 
// also non-embedded methods of the same order.
//
// Can be used to enable use of non-virtual calls for field, equation,
//   and stepper - potentially with inlined methods.
//
// Created: Josh Xie  June 2014 (supported by Google Summer of Code 2014 )
// Supervisors:  Sandro Wenzel, John Apostolakis (CERN)
//
// Adapted from G4CashKarpRKF45 class
// --------------------------------------------------------------------
// Original description (G4CashKarpRKF45):
// The Cash-Karp Runge-Kutta-Fehlberg 4/5 method is an embedded fourth
// order method (giving fifth-order accuracy) for the solution of an ODE.
// Two different fourth order estimates are calculated; their difference
// gives an error estimate. [ref. Numerical Recipes in C, 2nd Edition]
// Used to integrate the equations of motion of a particle in a field.
// Original Authors: J.Apostolakis, V.Grichine - 30.01.1997

#ifndef G4T_CASH_KARP_RKF45_HH
#define G4T_CASH_KARP_RKF45_HH

#include <cassert>

#include "G4LineSection.hh"
#include "G4MagIntegratorStepper.hh"

template <class T_Equation, unsigned int N = 6 >
class G4TCashKarpRKF45 : public G4MagIntegratorStepper
{
 public:

  G4TCashKarpRKF45(T_Equation* EqRhs, // G4int noIntegrationVariables = 6,
                    G4bool primary = true);

  virtual ~G4TCashKarpRKF45();

  inline void
  StepWithError(const G4double yInput[], // * __restrict__ yInput,
                const G4double dydx[],   // * __restrict__ dydx,
                G4double Step,
                G4double yOut[],         // * __restrict__ yOut,
                G4double yErr[] );       // * __restrict__ yErr);

  virtual void Stepper(const G4double yInput[],
                       const G4double dydx[],
                       G4double hstep,
                       G4double yOutput[],
                       G4double yError[]) override final;
  
  // __attribute__((always_inline))
  void RightHandSideInl( const G4double y[],  // * __restrict__  y,
                               G4double dydx[] ) // * __restrict__  dydx )
  {
    fEquation_Rhs->T_Equation::RightHandSide(y, dydx);
  }

  inline G4double DistChord() const override;

  inline G4int IntegratorOrder() const override { return 4; }

 private:
  G4TCashKarpRKF45(const G4TCashKarpRKF45&);
  G4TCashKarpRKF45& operator=(const G4TCashKarpRKF45&);
  // private copy constructor and assignment operator.

 private:
  G4double ak2[N], ak3[N], ak4[N], ak5[N], ak6[N], ak7[N], yTemp[N], yIn[N];
  // scratch space

  G4double fLastStepLength= 0.0;
  G4double* fLastInitialVector;
  G4double* fLastFinalVector;
  G4double* fLastDyDx;
  G4double* fMidVector;
  G4double* fMidError;
  // for DistChord calculations

  G4TCashKarpRKF45* fAuxStepper = nullptr;   
  // ... or G4TCashKarpRKF45<T_Equation, N>* fAuxStepper;
  T_Equation* fEquation_Rhs;
};

/////////////////////////////////////////////////////////////////////
//
// Constructor
//
template <class T_Equation, unsigned int N >
G4TCashKarpRKF45<T_Equation,N>::G4TCashKarpRKF45(T_Equation* EqRhs,
                                                 G4bool primary)
    : G4MagIntegratorStepper(dynamic_cast<G4EquationOfMotion*>(EqRhs), N )
    , fEquation_Rhs(EqRhs)
{
  if( dynamic_cast<G4EquationOfMotion*>(EqRhs) == nullptr )
  {
    G4Exception("G4TCashKarpRKF45: constructor", "GeomField0001",
                FatalException, "Equation is not an G4EquationOfMotion.");      
  }
    
  fLastInitialVector = new G4double[N];
  fLastFinalVector   = new G4double[N];
  fLastDyDx          = new G4double[N];
  
  fMidVector = new G4double[N];
  fMidError  = new G4double[N];
  
  if(primary)
  {
     fAuxStepper = new G4TCashKarpRKF45<T_Equation, N> (EqRhs, !primary);
  }
}

template <class T_Equation, unsigned int N >
G4TCashKarpRKF45<T_Equation,N>::~G4TCashKarpRKF45()
{
  delete[] fLastInitialVector;
  delete[] fLastFinalVector;
  delete[] fLastDyDx;
  delete[] fMidVector;
  delete[] fMidError;
  
  delete fAuxStepper;
}

//////////////////////////////////////////////////////////////////////
//
// Given values for n = 6 variables yIn[0,...,n-1]
// known  at x, use the fifth-order Cash-Karp Runge-
// Kutta-Fehlberg-4-5 method to advance the solution over an interval
// Step and return the incremented variables as yOut[0,...,n-1]. Also
// return an estimate of the local truncation error yErr[] using the
// embedded 4th-order method. The equation's method is called (inline)
// via RightHandSideInl(y,dydx), which returns derivatives dydx for y .
//
template <class T_Equation, unsigned int N >
inline void
G4TCashKarpRKF45<T_Equation,N>::StepWithError(const G4double* yInput,
                                              const G4double* dydx,
                                              G4double Step,
                                              G4double * yOut,
                                              G4double * yErr)
{
  // const G4double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875;

  const G4double b21 = 0.2, b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3,
                 b42 = -0.9, b43 = 1.2,
     
                 b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0,
                 b54 = 35.0 / 27.0,

                 b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0,
                 b63 = 575.0 / 13824.0, b64 = 44275.0 / 110592.0,
                 b65 = 253.0 / 4096.0,

                 c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 594.0,
                 c6 = 512.0 / 1771.0, dc5 = -277.0 / 14336.0;

  const G4double dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
                 dc4 = c4 - 13525.0 / 55296.0, dc6 = c6 - 0.25;

  // Initialise time to t0, needed when it is not updated by the integration.
  //       [ Note: Only for time dependent fields (usually electric)
  //                 is it neccessary to integrate the time.]
  // yOut[7] = yTemp[7]   = yIn[7];

  //  Saving yInput because yInput and yOut can be aliases for same array
  for(unsigned int i = 0; i < N; ++i)
  {
    yIn[i] = yInput[i];
  }
  // RightHandSideInl(yIn, dydx) ;              // 1st Step

  for(unsigned int i = 0; i < N; ++i)
  {
    yTemp[i] = yIn[i] + b21 * Step * dydx[i];
  }
  this->RightHandSideInl(yTemp, ak2);  // 2nd Step

  for(unsigned int i = 0; i < N; ++i)
  {
    yTemp[i] = yIn[i] + Step * (b31 * dydx[i] + b32 * ak2[i]);
  }
  this->RightHandSideInl(yTemp, ak3);  // 3rd Step
  
  for(unsigned int i = 0; i < N; ++i)
  {
    yTemp[i] = yIn[i] + Step * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
  }
  this->RightHandSideInl(yTemp, ak4);  // 4th Step
  
  for(unsigned int i = 0; i < N; ++i)
  {
    yTemp[i] = yIn[i] + Step * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] +
                                b54 * ak4[i]);
  }
  this->RightHandSideInl(yTemp, ak5);  // 5th Step
  
  for(unsigned int i = 0; i < N; ++i)
  {
    yTemp[i] = yIn[i] + Step * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] +
                                b64 * ak4[i] + b65 * ak5[i]);
  }
  this->RightHandSideInl(yTemp, ak6);  // 6th Step

  for(unsigned int i = 0; i < N; ++i)
  {
    // Accumulate increments with proper weights
     
    yOut[i] = yIn[i] +
       Step * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
  }
  for(unsigned int i = 0; i < N; ++i)
  {
    // Estimate error as difference between 4th and
    // 5th order methods
     
    yErr[i] = Step * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] +
                      dc5 * ak5[i] + dc6 * ak6[i]);
  }
  for(unsigned int i = 0; i < N; ++i)
  {
    // Store Input and Final values, for possible use in calculating chord
    fLastInitialVector[i] = yIn[i];
    fLastFinalVector[i]   = yOut[i];
    fLastDyDx[i]          = dydx[i];
  }
  // NormaliseTangentVector( yOut ); // Not wanted
  
  fLastStepLength = Step;
  
  return;
}

template <class T_Equation, unsigned int N >
inline G4double 
G4TCashKarpRKF45<T_Equation,N>::DistChord() const     
{
  G4double distLine, distChord;
  G4ThreeVector initialPoint, finalPoint, midPoint;
  
  // Store last initial and final points (they will be overwritten in
  // self-Stepper call!)
  initialPoint = G4ThreeVector(fLastInitialVector[0], fLastInitialVector[1],
                               fLastInitialVector[2]);
  finalPoint   = G4ThreeVector(fLastFinalVector[0], fLastFinalVector[1],
                               fLastFinalVector[2]);

  // Do half a step using StepNoErr
  
  fAuxStepper->G4TCashKarpRKF45::Stepper(fLastInitialVector, fLastDyDx,
                                         0.5 * fLastStepLength, fMidVector,
                                         fMidError);
  
  midPoint = G4ThreeVector(fMidVector[0], fMidVector[1], fMidVector[2]);
  
  // Use stored values of Initial and Endpoint + new Midpoint to evaluate
  //  distance of Chord
  
  if(initialPoint != finalPoint)
  {
    distLine  = G4LineSection::Distline(midPoint, initialPoint, finalPoint);
    distChord = distLine;
  }
  else
  {
    distChord = (midPoint - initialPoint).mag();
  }
  return distChord;
}

template <class T_Equation, unsigned int N >
inline void
G4TCashKarpRKF45<T_Equation,N>::Stepper(const G4double yInput[],
                                        const G4double dydx[],
                                        G4double Step,
                                        G4double yOutput[],
                                        G4double yError[])
{
  assert( yOutput != yInput );
  assert( yError  != yInput );
  
  StepWithError( yInput, dydx, Step, yOutput, yError);
}



#endif /* G4TCashKARP_RKF45_hh */
