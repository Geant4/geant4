// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidTorus.hh,v 1.2 1999-12-15 14:49:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __G4BREPSolidTorus
#define __G4BREPSolidTorus
#include "G4BREPSolid.hh"

class G4BREPSolidTorus: public G4BREPSolid
{
 public:
  G4BREPSolidTorus(const G4String ,
		   const G4ThreeVector&,
		   const G4ThreeVector&,
		   const G4ThreeVector&,
		   G4double,
		   G4double);
};

#endif
