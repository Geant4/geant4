// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ImplicitEuler.hh,v 1.2 1999-12-15 14:49:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// W. Wander <wwc@mit.edu> 12/09/97

#ifndef G4IMPLICITEULER_HH
#define G4IMPLICITEULER_HH
#include "G4MagErrorStepper.hh"

class G4ImplicitEuler: public G4MagErrorStepper
{

  public:
    G4ImplicitEuler(G4Mag_EqRhs *EqRhs, G4int numberOfVariables = 6) ;
   ~G4ImplicitEuler();

    void  DumbStepper(  const G4double y[] ,
		        const G4double dydx[] ,
		        const G4double h ,
			      G4double yout[] ) ;

    G4int     IntegratorOrder() { return 2 ; } ;

private: 

      G4int fNumberOfVariables ;
};

#endif /* G4IMPLICITEULER_HH */
