// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidCone.hh,v 1.4 2000-11-08 14:21:59 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BREPSolidCone
//
// Class description:
//
//  Definition of a generic BREP cone.
//
//  G4BREPSolidCone(const G4String& name,
//                  const G4ThreeVector& origin,
//                  const G4ThreeVector& axis,
//                  const G4ThreeVector& direction,
//                        G4double length,
//                        G4double radius,
//                        G4double large_radius)

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BREPSolidCone
#define __G4BREPSolidCone

#include "G4BREPSolid.hh"

class G4BREPSolidCone : public G4BREPSolid
{
 public:  // with description

  G4BREPSolidCone(const G4String& name, const G4ThreeVector& origin,
                  const G4ThreeVector& axis, const G4ThreeVector& direction,		   
                  G4double length, G4double radius, G4double large_radius);
    // Constructor

  ~G4BREPSolidCone();
    // Empty destructor.

  void Initialize();
    // Computes the bounding box for solids and surfaces.
    // Converts concave planes to convex.

  EInside Inside(register const G4ThreeVector& Pt) const;
    // Determines if the point Pt is inside, outside or on the surface
    // of the solid.

  G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const;
    // Returns the outwards pointing unit normal of the shape for the
    // surface closest to the point at offset p.

  G4double DistanceToIn(const G4ThreeVector& p) const;
    // Calculate the distance to the nearest surface of a shape from an
    // outside point p. The distance can be an underestimate.

  G4double DistanceToIn(register const G4ThreeVector& p, 
			register const G4ThreeVector& v) const;
    // Returns the distance along the normalised vector v to the shape,
    // from the point at offset p. If there is no intersection, returns
    // kInfinity. The first intersection resulting from `leaving' a
    // surface/volume is discarded. Hence, it is tolerant of points on
    // the surface of the shape.

  G4double DistanceToOut(register const G4ThreeVector& p, 
			 register const G4ThreeVector& v, 
			 const G4bool calcNorm=false, 
			 G4bool *validNorm=0, G4ThreeVector *n=0) const;
    // Returns the distance along the normalised vector v to the shape,
    // from a point at an offset p inside or on the surface of the shape.

  G4double DistanceToOut(const G4ThreeVector& p) const;
    // Calculates the distance to the nearest surface of a shape from an
    // inside point. The distance can be an underestimate.

 private:

  G4BREPSolidCone(const G4BREPSolidCone&);
  G4BREPSolidCone& operator=(const G4BREPSolidCone&);
    // Private copy constructor and assignment operator.

};

#endif
