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
// G4TMagErrorStepper
//
// Class description:
//
// Templated version of G4MagErrorStepper
//
//
// Created: Josh Xie  (supported by Google Summer of Code 2014 )
// Supervisors:  Sandro Wenzel, John Apostolakis (CERN)
// Adapted from G4G4TMagErrorStepper class
// --------------------------------------------------------------------
#ifndef G4TMAG_ERROR_STEPPER_HH
#define G4TMAG_ERROR_STEPPER_HH

#include "G4Types.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ThreeVector.hh"
#include "G4LineSection.hh"

template <class T_Stepper, class T_Equation, unsigned int N>
class G4TMagErrorStepper : public G4MagIntegratorStepper
{
 public:  // with description
  G4TMagErrorStepper(T_Equation* EqRhs, G4int numberOfVariables,
                     G4int numStateVariables = 12)
    : G4MagIntegratorStepper(EqRhs, numberOfVariables, numStateVariables)
    , fEquation_Rhs(EqRhs)
  {
    // G4int nvar = std::max(this->GetNumberOfVariables(), 8);
  }

  virtual ~G4TMagErrorStepper() { ; }

  inline void RightHandSide(G4double y[], G4double dydx[])
  {
    fEquation_Rhs->T_Equation::RightHandSide(y, dydx);
  }

  inline void Stepper(const G4double yInput[], const G4double dydx[],
                      G4double hstep, G4double yOutput[], G4double yError[]) override final;
 
  inline G4double DistChord() const override final;

 private:
  G4TMagErrorStepper(const G4TMagErrorStepper&);
  G4TMagErrorStepper& operator=(const G4TMagErrorStepper&);
  // Private copy constructor and assignment operator.

 private:
  // STATE
  G4ThreeVector fInitialPoint, fMidPoint, fFinalPoint;
  // Data stored in order to find the chord

  // Dependent Objects, owned --- part of the STATE
  G4double yInitial[N < 8 ? 8 : N];
  G4double yMiddle[N < 8 ? 8 : N];
  G4double dydxMid[N < 8 ? 8 : N];
  G4double yOneStep[N < 8 ? 8 : N];
  // The following arrays are used only for temporary storage
  // they are allocated at the class level only for efficiency -
  // so that calls to new and delete are not made in Stepper().

  T_Equation* fEquation_Rhs;
};

// ------------   Implementation -----------------------

template <class T_Stepper, class T_Equation, unsigned int N >
void G4TMagErrorStepper<T_Stepper,T_Equation,N>::
Stepper(const G4double yInput[],
        const G4double dydx[],
              G4double hstep,
              G4double yOutput[],
              G4double yError[])
// The stepper for the Runge Kutta integration. The stepsize
// is fixed, with the Step size given by hstep.
// Integrates ODE starting values y[0 to N].
// Outputs yout[] and its estimated error yerr[].
{
  const unsigned int maxvar = GetNumberOfStateVariables();
  
  //  Saving yInput because yInput and yOutput can be aliases for same array
  for(unsigned int i = 0; i < N; ++i)
     yInitial[i] = yInput[i];
  yInitial[7] =
     yInput[7];  // Copy the time in case ... even if not really needed
  yMiddle[7]  = yInput[7];  // Copy the time from initial value
  yOneStep[7] = yInput[7];  // As it contributes to final value of yOutput ?
  // yOutput[7] = yInput[7];  // -> dumb stepper does it too for RK4
  for(unsigned int i = N; i < maxvar; ++i)
     yOutput[i] = yInput[i];

  G4double halfStep = hstep * 0.5;
  
  // Do two half steps
  
  static_cast<T_Stepper*>(this)->DumbStepper(yInitial, dydx, halfStep,
                                               yMiddle);
  this->RightHandSide(yMiddle, dydxMid);
  static_cast<T_Stepper*>(this)->DumbStepper(yMiddle, dydxMid, halfStep,
                                             yOutput);
  
  // Store midpoint, chord calculation
  
  fMidPoint = G4ThreeVector(yMiddle[0], yMiddle[1], yMiddle[2]);
  
  // Do a full Step
  static_cast<T_Stepper*>(this)->DumbStepper(yInitial, dydx, hstep, yOneStep);
  for(unsigned int i = 0; i < N; ++i)
  {
    yError[i] = yOutput[i] - yOneStep[i];
    yOutput[i] +=
       yError[i] *
       T_Stepper::IntegratorCorrection;  // Provides accuracy increased
    // by 1 order via the
    // Richardson Extrapolation
  }

  fInitialPoint = G4ThreeVector(yInitial[0], yInitial[1], yInitial[2]);
  fFinalPoint   = G4ThreeVector(yOutput[0], yOutput[1], yOutput[2]);
  
  return;
}


template <class T_Stepper, class T_Equation, unsigned int N >
inline G4double
G4TMagErrorStepper<T_Stepper,T_Equation,N>::DistChord() const
{
  // Estimate the maximum distance from the curve to the chord
  //
  //  We estimate this using the distance of the midpoint to
  //  chord (the line between
  //
  //  Method below is good only for angle deviations < 2 pi,
  //   This restriction should not a problem for the Runge cutta methods,
  //   which generally cannot integrate accurately for large angle deviations.
  G4double distLine, distChord;

  if(fInitialPoint != fFinalPoint)
  {
    distLine = G4LineSection::Distline(fMidPoint, fInitialPoint, fFinalPoint);
    // This is a class method that gives distance of Mid
    //  from the Chord between the Initial and Final points.
    
    distChord = distLine;
  }
  else
  {
     distChord = (fMidPoint - fInitialPoint).mag();
  }
  
  return distChord;
}

#endif /* G4TMagErrorStepper_HH */
