//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4MagIntegratorStepper.hh,v 1.9 2002-11-20 18:09:22 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// - 20.11.02  J. Apostolakis: Added new "State" elements

#ifndef G4MAGIntegratorSTEPPER
#define G4MAGIntegratorSTEPPER

#include "globals.hh"
#include "G4EquationOfMotion.hh"

class G4MagIntegratorStepper
{
  public:  // with description

     G4MagIntegratorStepper(G4EquationOfMotion *Equation, 
			    G4int              numIntegrationVariables,
			    G4int              numStateVariables=12);
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

     inline void NormaliseTangentVector( G4double vec[6] );
       // Simple utility function to (re)normalise 'unit velocity' vector.

     inline void RightHandSide( const double y[], double dydx[] );   
       // Utility method to supply the standard Evaluation of the
       // Right Hand side of the associated equation.

     inline G4int  GetNumberOfVariables() const;
       // Get the number of variables that the stepper will integrate over.

     // void   SetNumberOfVariables(G4int newNo);  // Dangerous & obsolete ...

     inline G4int  GetNumberOfStateVariables() const;
       // Get the number of variables of state variables (>= above, integration)

     virtual G4int IntegratorOrder() const = 0;
       // Returns the order of the integrator
       // i.e. its error behaviour is of the order O(h^order).

     inline G4EquationOfMotion *GetEquationOfMotion(); 
       // As some steppers (eg RKG3) require other methods of Eq_Rhs
       // this function allows for access to them.

  public:  // without description

#if 0
     void
     SetChargeAndMomentum( G4double particleCharge, // in e+ units
			   G4double MomentumXc)
       //  Supply the standard Evaluation of the Right Hand side 
       //   of the associated equation.
       {theEquation_Rhs -> SetChargeAndMomentum(particleCharge, MomentumXc);}
#endif 

  private:
  
     G4MagIntegratorStepper(const G4MagIntegratorStepper&);
     G4MagIntegratorStepper& operator=(const G4MagIntegratorStepper&);
       // Private copy constructor and assignment operator.

  private:

     G4EquationOfMotion *fEquation_Rhs;
     const G4int  fNoIntegrationVariables;  // Number of Variables in integration
     const G4int  fNoStateVariables;        // Number required for FieldTrack
     // const G4int  fNumberOfVariables;
};

#include  "G4MagIntegratorStepper.icc"

#endif  /* G4MAGIntegratorSTEPPER */
