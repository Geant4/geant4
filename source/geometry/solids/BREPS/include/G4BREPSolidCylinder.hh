// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidCylinder.hh,v 2.1 1998/10/20 16:31:07 broglia Exp $
// GEANT4 tag $Name: geant4-00 $
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
