// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SimpleRunge.hh,v 1.3 2000-04-27 09:14:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4SimpleRunge
//
// Class description:
//
// Simple Runge:
//
//        x_1 = x_0 + h * ( dx( t_0+h/2, x_0 + h/2 * dx( t_0, x_0) ) )
//
// second order solver.
// Take the derivative at a position to be assumed at the middle of the
// Step and add it to the current position.

// History:
// - Created. W.Wander <wwc@mit.edu>, 12/09/97

#ifndef G4SIMPLERUNGE_HH
#define G4SIMPLERUNGE_HH

#include "G4MagErrorStepper.hh"

class G4SimpleRunge : public G4MagErrorStepper
{

  public:  // with description

    G4SimpleRunge(G4Mag_EqRhs *EqRhs, G4int numberOfVariables = 6) ;
   ~G4SimpleRunge();
      // Constructor and destructor.

    void DumbStepper( const G4double y[],
		      const G4double dydx[],
		      const G4double h,
			    G4double yout[]);

  public:  // without description
  
    G4int IntegratorOrder() { return 2; }

  private:

    G4int fNumberOfVariables ;

    G4double* dydxTemp;
    G4double* dydxTemp2;
    G4double* yTemp;
      // scratch space    
};

#endif /* G4SIMPLERUNGE_HH */
