// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidBox.hh,v 1.2 1999-12-15 14:49:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __G4BREPSolidBOX
#define __G4BREPSolidBOX
#include "G4BREPSolid.hh"
#include "G4RotationMatrix.hh"

class G4BREPSolidBox: public G4BREPSolid
{
public:

  G4BREPSolidBox(G4String,const G4Point3D&, const G4Point3D&, 
		 const G4Point3D&, const G4Point3D&, const G4Point3D&, 
		 const G4Point3D&, const G4Point3D&, const G4Point3D& );       

  EInside Inside(register const G4ThreeVector&) const;

private:

  G4RotationMatrix Rotation;

};

#endif
