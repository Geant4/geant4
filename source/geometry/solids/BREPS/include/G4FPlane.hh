// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FPlane.hh,v 1.8 2000-08-28 08:57:46 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4FPlane
//
// Class Description:
//   
// A G4FPlane is a plane created by 3 points or by an origin, an axis and 
// a direction. The plane created is a G4Plane, where his coefficient a, b, 
// c and d are stored. The equation of the plane is:
//                ax + by + cz = d
// 
// This class contain 2 intersection functions :
//       - closest intersection 
//       - intersection by a ray

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, S.Giani, G.Cosmo.
// ----------------------------------------------------------------------
//
// History
// -------
// - SurfaceNormal always returns the direction of NormalX, always containing
//   the correct orientation for all faces (S.Giani).
// - Addition of default argument sense = 1 in the second constructor (S.Giani).
// - The constructor using iVec now properly stores both the internal and
//   external boundaries in the bounds vector (S.Giani).
// - Proper initialization of sameSense in both the constructors (S.Giani). 
// - Addition of third argument (sense) in the second constructor to ensure
//   consistent setting of the normal in all the client code (S.Giani).
// - Proper use of the tolerance in the Intersect function (S.Giani).
// ----------------------------------------------------------------------
#ifndef __PLANESURFACE_H
#define __PLANESURFACE_H

#include "G4Axis2Placement3D.hh"
#include "G4Plane.hh"
#include "G4Surface.hh"


class G4FPlane : public G4Surface
{

public:  // with description

  G4FPlane();
  ~G4FPlane();
    // Default constructor & destructor.

  G4FPlane( const G4Vector3D& direction, 
	    const G4Vector3D& axis     ,	   
	    const G4Point3D&  Pt0       );
    // Normal constructor.

  G4FPlane(const G4Point3DVector* pVec, 
	   const G4Point3DVector* iVec= 0,
	   G4int sense = 1);
    // Constructor used by G4BREPSolidBox and G4BREPSolidPolyhedra.

  G4int Intersect(const G4Ray& G4Rayref);
    // Calculates the intersection of the plane and a ray.

  void CalcBBox();
    // Calculates bounding box.

  void Project();
    // Computes the projection of the plane.
  
  inline G4int GetConvex();
    // Return plane's convexity, if so.

  inline G4int GetNumberOfPoints();
    // Gets the number of the points on the surface boundary.

  inline G4Point3D GetSrfPoint();
    // Gets the location point on the surface.

  inline const G4Point3D& GetPoint(G4int Count);
    // Gets a surface boundary point.

  void CalcNormal();
    // Computes normal to surface.

  inline G4Vector3D SurfaceNormal(const G4Point3D& Pt) const;
    // Returns normal to surface.

  inline const char* Name() const;
    // Returns the type identifier.

  G4double ClosestDistanceToPoint(const G4Point3D& Pt);
    // Returns the closest distance from point Pt.

  G4double HowNear( const G4Vector3D& x ) const ;
    // Computes the shortest distance from the point x to the G4FPlane.
    // The distance will always be positive.

  inline G4Axis2Placement3D GetPplace() const;
  inline G4Plane GetPplane() const;
    // Accessors to geometrical data.

public:  // without description

  inline G4int MyType() const;
    // Returns the shape type (used in G4BREPSolid).
  
  G4int IsConvex(); 
    // Returns -1.  (?)

  inline void Deactivate();
    // Deactive, used in G4Surface.

  inline G4Ray* Norm();
    // Returns the normal (used in BREPSolid).

  G4Point3D hitpoint;
    // Hit point of the ray on the surface.

protected:

  void InitBounded();
  
private:

  inline G4int Sign(G4double a);

private:

  G4Axis2Placement3D pplace;
  G4Plane            Pl;
  G4Ray              *NormalX;
  G4int              Convex;
  G4SurfaceBoundary* projectedBoundary;

};

#include "G4FPlane.icc"

#endif
