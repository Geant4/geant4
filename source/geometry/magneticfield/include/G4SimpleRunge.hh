// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SimpleRunge.hh,v 1.1 1999-01-07 16:07:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// W. Wander <wwc@mit.edu> 12/09/97

#ifndef G4SIMPLERUNGE_HH
#define G4SIMPLERUNGE_HH
#include "G4MagErrorStepper.hh"

class G4SimpleRunge: public G4MagErrorStepper
{

  public:
    G4SimpleRunge(G4Mag_EqRhs *EqRhs, G4int numberOfVariables = 6) ;
   ~G4SimpleRunge();

    void  DumbStepper(  const G4double y[],
		    const G4double dydx[],
		    const G4double h,
			  G4double yout[]);

    G4int     IntegratorOrder() { return 2; };

private:

    G4int fNumberOfVariables ;

    // scratch space    
    G4double* dydxTemp;
    G4double* dydxTemp2;
    G4double* yTemp;
};

#endif /* G4SIMPLERUNGE_HH */
