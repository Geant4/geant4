// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidTorus.hh,v 1.3 2000-08-28 08:57:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BREPSolidTorus
//
// Class description:
// 
// Definition of a generic BREP torus.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BREPSolidTorus
#define __G4BREPSolidTorus

#include "G4BREPSolid.hh"

class G4BREPSolidTorus : public G4BREPSolid
{

 public:  // with description

  G4BREPSolidTorus(const G4String& name,
		   const G4ThreeVector& origin,
		   const G4ThreeVector& axis,
		   const G4ThreeVector& direction,
		   G4double MinRadius,
		   G4double MaxRadius);
    // Constructor defining the torus through its surfaces.
};

#endif
