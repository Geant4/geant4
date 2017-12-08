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
// $Id: G4MagIntegratorStepper.hh 107059 2017-11-01 14:58:16Z gcosmo $
//
//
// class G4MagIntegratorStepper
//
// Class description:
//
// Abstract base class for integrator of particle's equation of motion,
// used in tracking in space dependent magnetic field
//
//  A Stepper must integrate over                NumberOfVariables elements,
//   and also copy (from input to output) any of NoStateVariables  
//   not included in the NumberOfVariables.  
// 
//  So it is expected that NoStateVariables >= NumberOfVariables

// History:
// - 15.01.97  J. Apostolakis (J.Apostolakis@cern.ch)
// --------------------------------------------------------------------

#ifndef G4MAGIntegratorSTEPPER
#define G4MAGIntegratorSTEPPER

#include "G4Types.hh"
#include "G4EquationOfMotion.hh"

class G4MagIntegratorStepper
{
  public:  // with description

     G4MagIntegratorStepper(G4EquationOfMotion *Equation, 
                            G4int              numIntegrationVariables,
                            G4int              numStateVariables=12,
                            bool               isFSAL= false
                            // , G4int       methodOrder
                           );
     virtual ~G4MagIntegratorStepper();
       // Constructor and destructor. No actions.

     virtual  void  Stepper(  const G4double y[],
                              const G4double dydx[],
                                    G4double h,
                                    G4double yout[],
                                    G4double yerr[]  ) = 0 ;
       // The stepper for the Runge Kutta integration.
       // The stepsize is fixed, with the Step size given by h.
       // Integrates ODE starting values y[0 to 6].
       // Outputs yout[] and its estimated error yerr[].

     virtual  G4double  DistChord() const = 0;
       // Estimate the maximum distance of a chord from the true path
       // over the segment last integrated.

     virtual void ComputeRightHandSide( const G4double y[], G4double dydx[] ); 
       // Must compute the RightHandSide as in the method below
       // Optionally can cache the input y[] and the dydx[] values computed.

     inline void NormaliseTangentVector( G4double vec[6] );
       // Simple utility function to (re)normalise 'unit velocity' vector.

     inline void NormalisePolarizationVector( G4double vec[12] );
       // Simple utility function to (re)normalise 'unit spin' vector.

     inline void RightHandSide( const double y[], double dydx[] ) const;
       // Utility method to supply the standard Evaluation of the
       // Right Hand side of the associated equation.

     inline G4int  GetNumberOfVariables() const;
       // Get the number of variables that the stepper will integrate over.

     inline G4int  GetNumberOfStateVariables() const;
       // Get the number of variables of state variables (>= above, integration)

     virtual G4int IntegratorOrder() const = 0;
       // Returns the order of the integrator
       // i.e. its error behaviour is of the order O(h^order).

     G4int IntegrationOrder() { return fIntegrationOrder; }
       //  Replacement method - using new data member
   
     inline G4EquationOfMotion *GetEquationOfMotion(); 
       // As some steppers (eg RKG3) require other methods of Eq_Rhs
       // this function allows for access to them.
     inline void SetEquationOfMotion(G4EquationOfMotion* newEquation); 

     inline unsigned long GetfNoRHSCalls(){ return fNoRHSCalls; }
     // void IncrementRHSCalls() { fNoRHSCalls++; }
     inline void ResetfNORHSCalls(){ fNoRHSCalls = 0; }
       // Count number of calls to RHS method(s)

     bool IsFSAL() const { return fIsFSAL; }
   
  protected:
     void SetIntegrationOrder(int order) { fIntegrationOrder= order; }
     void SetFSAL( bool flag= true) { fIsFSAL= flag; }
   
  private:
  
     G4MagIntegratorStepper(const G4MagIntegratorStepper&);
     G4MagIntegratorStepper& operator=(const G4MagIntegratorStepper&);
       // Private copy constructor and assignment operator.

  private:

     G4EquationOfMotion *fEquation_Rhs;
     const G4int  fNoIntegrationVariables;  // Number of Variables in integration
     const G4int  fNoStateVariables;        // Number required for FieldTrack
     // const G4int  fNumberOfVariables;

     // Counter for calls to RHS method
     mutable unsigned long fNoRHSCalls;

     // Parameters of a RK method -- must be shared by all steppers of a type
     // -- Invariants for a class
     /* const */ int     fIntegrationOrder;  // all ClassicalRK4 steppers are 4th order
     /* const */ bool    fIsFSAL;            // Depends on RK method & implementation
};

#include  "G4MagIntegratorStepper.icc"

#endif  /* G4MAGIntegratorSTEPPER */
