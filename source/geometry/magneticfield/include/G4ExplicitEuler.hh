// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExplicitEuler.hh,v 1.1 1999-01-07 16:07:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// W. Wander <wwc@mit.edu> 12/09/97

#ifndef G4EXPLICITEULER_HH
#define G4EXPLICITEULER_HH
#include "G4MagErrorStepper.hh"

class G4ExplicitEuler: public G4MagErrorStepper
{

  public:
    G4ExplicitEuler(G4Mag_EqRhs *EqRhs, G4int numberOfVariables = 6) ;
   ~G4ExplicitEuler();

    void  DumbStepper(  const G4double y[],
		        const G4double dydx[],
		        const G4double h,
			      G4double yout[]);

    G4int     IntegratorOrder() { return 1; };

private:

    G4int fNumberOfVariables ;
};

#endif /* G4EXPLICITEULER_HH */
