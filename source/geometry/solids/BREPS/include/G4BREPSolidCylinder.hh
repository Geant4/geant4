// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidCylinder.hh,v 1.3 2000-08-28 08:57:41 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BREPSolidCylinder
//
// Class description:
//
// Definition of a generic BREP cylinder.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BREPSolidCylinder
#define __G4BREPSolidCylinder

#include "G4BREPSolid.hh"

class G4BREPSolidCylinder : public G4BREPSolid
{

public: // with description

  G4BREPSolidCylinder(const G4String& name,
		      const G4ThreeVector& origin,
		      const G4ThreeVector& axis,
		      const G4ThreeVector& direction,		      
		      G4double radius,
		      G4double length);
    // Constructor defining the cylinder through planes
    // and circular surfaces.
};

#endif
