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
// G4MagIntegratorStepper
//
// Class description:
//
// Abstract base class for integrator of particle's equation of motion,
// used in tracking in space dependent magnetic field
//
// A Stepper must integrate over NumberOfVariables elements,
// and also copy (from input to output) any of NoStateVariables  
// not included in the NumberOfVariables.  
// 
// So it is expected that NoStateVariables >= NumberOfVariables

// Author: J.Apostolakis, CERN - 15.01.1997
// --------------------------------------------------------------------
#ifndef G4MAGINTEGRATORSTEPPER_HH
#define G4MAGINTEGRATORSTEPPER_HH

#include "G4Types.hh"
#include "G4EquationOfMotion.hh"
#include "G4VIntegrationDriver.hh"
#include "G4IntegrationDriver.hh"

class G4VIntegrationDriver;

class G4MagIntegratorStepper
{
  public:  // with description

     G4MagIntegratorStepper(G4EquationOfMotion* Equation, 
                            G4int               numIntegrationVariables,
                            G4int               numStateVariables = 12,
                            G4bool              isFSAL = false );

     virtual ~G4MagIntegratorStepper() = default;
       // Constructor and destructor. No actions.

     G4MagIntegratorStepper(const G4MagIntegratorStepper&) = delete;
     G4MagIntegratorStepper& operator=(const G4MagIntegratorStepper&) = delete;

     virtual void Stepper( const G4double y[],
                           const G4double dydx[],
                                 G4double h,
                                 G4double yout[],
                                 G4double yerr[] ) = 0;
       // The stepper for the Runge Kutta integration.
       // The stepsize is fixed, with the Step size given by h.
       // Integrates ODE starting values y[0 to 6].
       // Outputs yout[] and its estimated error yerr[].

     virtual G4double DistChord() const = 0;
       // Estimate the maximum distance of a chord from the true path
       // over the segment last integrated.

     inline void NormaliseTangentVector( G4double vec[6] );
       // Simple utility function to (re)normalise 'unit velocity' vector.

     inline void NormalisePolarizationVector( G4double vec[12] );
       // Simple utility function to (re)normalise 'unit spin' vector.

     inline void RightHandSide( const G4double y[], G4double dydx[] ) const;
       // Utility method to supply the standard Evaluation of the
       // Right Hand side of the associated equation.

     inline void RightHandSide( const G4double y[],
                                      G4double dydx[],
                                      G4double field[] ) const;
       // Calculate dydx and field at point y. 

     inline G4int  GetNumberOfVariables() const;
       // Get the number of variables that the stepper will integrate over.

     inline G4int  GetNumberOfStateVariables() const;
       // Get the number of variables of state variables (>= above, integration)

     virtual G4int IntegratorOrder() const = 0;
       // Returns the order of the integrator
       // i.e. its error behaviour is of the order O(h^order).

     inline G4int IntegrationOrder();
       // Replacement method - using new data member
   
     inline G4EquationOfMotion* GetEquationOfMotion();
     inline const G4EquationOfMotion* GetEquationOfMotion() const;
       // As some steppers (eg RKG3) require other methods of Eq_Rhs
       // this function allows for access to them.

     inline void SetEquationOfMotion(G4EquationOfMotion* newEquation); 

     inline unsigned long GetfNoRHSCalls();
     inline void ResetfNORHSCalls();
       // Count number of calls to RHS method(s)

     inline G4bool IsFSAL() const;

     // TODO - QSS
     inline G4bool isQSS() const { return fIsQSS; }
     void   SetIsQSS(G4bool val){ fIsQSS= val;}

protected:

     inline void SetIntegrationOrder(G4int order);
     inline void SetFSAL(G4bool flag = true);
   
  private:

     G4EquationOfMotion* fEquation_Rhs = nullptr;
     const G4int fNoIntegrationVariables = 0; // Variables in integration
     const G4int fNoStateVariables = 0;       // Number required for FieldTrack

     mutable unsigned long fNoRHSCalls = 0UL;
       // Counter for calls to RHS method

     // Parameters of a RK method -- must be shared by all steppers of a type
     // -- Invariants for a class
 
     G4int fIntegrationOrder = -1;  // must be set by stepper !!!
      // All ClassicalRK4 steppers are 4th order
     G4bool fIsFSAL = false;
      // Depends on RK method & implementation
     G4bool fIsQSS = false;

};

#include  "G4MagIntegratorStepper.icc"

#endif  /* G4MAGIntegratorSTEPPER */
