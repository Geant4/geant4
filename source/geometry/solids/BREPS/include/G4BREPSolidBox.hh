// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidBox.hh,v 1.3 2000-08-28 08:57:41 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BREPSolidBox
//
// Class description:
//
// Definition of a generic BREP box.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BREPSolidBOX
#define __G4BREPSolidBOX

#include "G4BREPSolid.hh"
#include "G4RotationMatrix.hh"

class G4BREPSolidBox : public G4BREPSolid
{

public: // with description

  G4BREPSolidBox(const G4String&, const G4Point3D&, const G4Point3D&, 
		 const G4Point3D&, const G4Point3D&, const G4Point3D&, 
		 const G4Point3D&, const G4Point3D&, const G4Point3D& );       
    // Constructor defining the box through its planes.

  EInside Inside(register const G4ThreeVector& Pt) const;
    // Determines if the point Pt is inside, outside or on the surface
    // of the solid.

private:

  G4RotationMatrix Rotation;

};

#endif
