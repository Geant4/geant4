// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidSphere.hh,v 1.1 1999-01-07 16:07:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __G4BREPSolidSphere
#define __G4BREPSolidSphere
#include "G4BREPSolid.hh"

class G4BREPSolidSphere: public G4BREPSolid
{
 public:
  G4BREPSolidSphere(const G4String,
		    const G4Vector3D&,
		    const G4Vector3D&, 
		    const G4Vector3D&,
		    G4double);

  inline void SphReset()const
  {
    ((G4BREPSolidSphere*)this)->active=1;
  }

  EInside Inside(register const G4ThreeVector&) const;
  G4ThreeVector SurfaceNormal(const G4ThreeVector&) const;

  G4double DistanceToIn(const G4ThreeVector&) const;
  G4double DistanceToIn(register const G4ThreeVector&, 
			register const G4ThreeVector&) const;

  G4double DistanceToOut(register const G4ThreeVector&, 
			 register const G4ThreeVector&, 
			 const G4bool calcNorm=false, 
			 G4bool *validNorm=0, G4ThreeVector *n=0) const;
  G4double DistanceToOut(const G4ThreeVector&) const;
};

#endif






