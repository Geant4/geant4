// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HelixImplicitEuler.hh,v 1.1 1999-01-07 16:07:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// W. Wander <wwc@mit.edu> 03/11/98

#ifndef G4HELIXIMPLICITEULER_HH
#define G4HELIXIMPLICITEULER_HH
#include "G4MagHelicalStepper.hh"

class G4HelixImplicitEuler: public G4MagHelicalStepper
{

public:
  G4HelixImplicitEuler(G4Mag_EqRhs *EqRhs): G4MagHelicalStepper(EqRhs){};
  ~G4HelixImplicitEuler(){};
  
  void  DumbStepper(  const G4double y[],
		      const G4double dydx[],
		      const G4double h,
		      G4double yout[]);
  
  G4int     IntegratorOrder() { return 2; };
};

#endif /* G4HELIXIMPLICITEULER_HH */
