// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidTorus.cc,v 1.4 2000-11-08 14:22:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BREPSolidTorus.cc
//
// ----------------------------------------------------------------------

#include "G4BREPSolidTorus.hh"
#include "G4ToroidalSurface.hh"

G4BREPSolidTorus::G4BREPSolidTorus(const    G4String& name,
				   const    G4ThreeVector& origin,
				   const    G4ThreeVector& axis,
				   const    G4ThreeVector& direction,
				   G4double MinRadius,
				   G4double MaxRadius)
  : G4BREPSolid(name)
{
  SurfaceVec = new G4Surface*[1];
  SurfaceVec[0] = new G4ToroidalSurface( origin, axis, direction,
					 MinRadius, MaxRadius);
  nb_of_surfaces = 1;
  active = 1;
  Initialize();
}

G4BREPSolidTorus::~G4BREPSolidTorus()
{
}
