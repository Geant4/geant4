// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MagIntegratorStepper.hh,v 1.2 1999-02-12 12:28:40 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Abstract base class (ie Interface)
// -------------------
//    for integrator of particle's equation of motion,
//    used in tracking in space dependent magnetic field
// -----------------------------------------------------
//
// History:
// 15.01.97  J.Apostolakis (J.Apostolakis@cern.ch)

#ifndef G4MAGIntegratorSTEPPER
#define G4MAGIntegratorSTEPPER
#include "globals.hh"
#include "G4Mag_EqRhs.hh"

class G4MagIntegratorStepper
{
  public:

     G4MagIntegratorStepper(G4Mag_EqRhs *EqRhs, G4int num_variables);
     ~G4MagIntegratorStepper(){} ;

     //  "Key" methods
     // ---------------
     //   The stepper for the Runge Kutta integration. The stepsize 
     // is fixed, with the Step size given by h.
     //  Integrates ODE starting values y[0 to 6 ]
     // Outputs yout[] and its estimated error yerr[].

     virtual  void  Stepper(  const G4double y[],
			      const G4double dydx[],
			      const G4double h,
				    G4double yout[],
				    G4double yerr[]  ) = 0 ;

     //  Estimate the maximum distance of chord from true path over
     //    segment last integrated.

     virtual  G4double  DistChord() const = 0; 


     // Utility methods
     // ---------------
     //  Simple function to (re)normalise 'unit velocity' vector
     //
     void NormaliseTangentVector( G4double vec[6] );

     //  Supply the standard Evaluation of the Right Hand side 
     //   of the associated equation.
     //
     virtual void RightHandSide( const  double y[], double dydx[] );   
                                        // FIXME : not virtual   JA 10/2/99

     //  Get/Set the number of variables that the stepper will compile over
     G4int  GetNumberOfVariables();
     void   SetNumberOfVariables(G4int newNo);

#if 0
     //  Supply the standard Evaluation of the Right Hand side 
     //   of the associated equation.
     void
     SetChargeAndMomentum( const G4double particleCharge, // in e+ units
			   const G4double MomentumXc)
       { theEquation_Rhs -> SetChargeAndMomentum(particleCharge, MomentumXc);}
#endif 

     // returns the order of the integrator
     // i.e. its error behaviour is of the order O(h^order)
     virtual G4int IntegratorOrder() = 0;
  
     //  As some steppers (eg RKG3) require other methods of Eq_Rhs
     //    the next function allows for access to them.
     G4EquationOfMotion *GetEquationOfMotion() const; 

  private:

     G4EquationOfMotion *fEquation_Rhs;
     G4int              fNumberOfVariables;

};

#include  "G4MagIntegratorStepper.icc"

#endif  /* G4MAGIntegratorSTEPPER */
