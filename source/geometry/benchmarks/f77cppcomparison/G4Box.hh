// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Box.hh,v 1.2 1999-12-15 14:49:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4Box
//
// A Box is a cuboid of given half lengths dx,dy,dz. The Box is
// centred on the origin with sides parallel to the x/y/z axes.
//
// Member functions:
//
// As inherited from G4VSolid +
//
// G4Box(const G4String& pName,const G4double pX,
//       const G4double pY,const G4double pZ)
//   Construct a box with name, and half lengths pX,pY,pZ
//
// G4double GetXHalfLength() const
// G4double GetYHalfLength() const
// G4double GetZHalfLength() const
//
//   Return the respective parameter
//
// Protected:
// 
// G4ThreeVectorList*
// CreateRotatedVertices(const G4Transform& pTransform) const
// 
//   Create the List of transformed vertices in the format required
//   for G4VSolid:: ClipCrossSection and ClipBetweenSections.
//
// Member Data:
//
// fDx,fDy,fDz - The box's half-widths
//
// History:
// 30.06.95 P.Kent Converted from source code developed end 94
// 18.07.95 J.Allison Added virtual function Wireframe.
// 31.07.95 J.Allison Added virtual function DispatchWireframe.

#ifndef G4BOX_HH
#define G4BOX_HH

#include "G4ThreeVector.hh"

class G4Box
{
public:
    G4Box(const G4double pX,
	  const G4double pY,const G4double pZ);

// Access functions
    G4double GetXHalfLength() const
    {
	return fDx;
    }
    
    G4double GetYHalfLength() const
    {
	return fDy;
    }
    
    G4double GetZHalfLength() const
    {
	return fDz;
    }
    
    EInside Inside(const G4ThreeVector& p) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

    G4double DistanceToIn(const G4ThreeVector& p,const G4ThreeVector& v) const;
    G4double DistanceToIn(const G4ThreeVector& p) const;
    G4double DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
			   const G4bool calcNorm=false,
			   G4bool *validNorm=0,G4ThreeVector *n=0) const;
    G4double DistanceToOut(const G4ThreeVector& p) const;

protected:
    G4double fDx,fDy,fDz;
};

#endif




