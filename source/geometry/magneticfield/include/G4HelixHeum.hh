// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HelixHeum.hh,v 1.1.10.1 1999/12/07 20:48:00 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
//
// W. Wander <wwc@mit.edu> 03/11/98

#ifndef G4HELIXHEUM_HH
#define G4HELIXHEUM_HH
#include "G4MagHelicalStepper.hh"

class G4HelixHeum: public G4MagHelicalStepper
{

public:
  G4HelixHeum(G4Mag_EqRhs *EqRhs): G4MagHelicalStepper(EqRhs){};
  ~G4HelixHeum(){};
  
  void  DumbStepper(  const G4double y[],
		      const G4double dydx[],
		      const G4double h,
		      G4double yout[]);
  
  G4int     IntegratorOrder() { return 2; };
};

#endif /* G4HELIXHEUM_HH */
