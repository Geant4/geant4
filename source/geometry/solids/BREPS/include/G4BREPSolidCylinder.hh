// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidCylinder.hh,v 1.1.10.1 1999/12/07 20:48:16 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
#ifndef __G4BREPSolidCylinder
#define __G4BREPSolidCylinder
#include "G4BREPSolid.hh"

class G4BREPSolidCylinder : public G4BREPSolid
{
public:
  G4BREPSolidCylinder(G4String name,
		      const G4ThreeVector&,
		      const G4ThreeVector&,
		      const G4ThreeVector&,		      
		      const G4double&,
		      const G4double&);
};

#endif
