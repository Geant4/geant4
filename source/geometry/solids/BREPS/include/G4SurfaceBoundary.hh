// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SurfaceBoundary.hh,v 1.4 2000-08-28 15:00:34 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4SurfaceBoundary
//
// Class description:
// 
// Definition of a surface boundary.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef included_G4SurfaceBoundary
#define included_G4SurfaceBoundary

#include "G4Point3D.hh"
#include "G4Point3DVector.hh"
#include "G4Vector3D.hh"
#include "G4Transform3D.hh"
#include "G4Curve.hh"
#include "G4CurveVector.hh"
#include "G4CurveRayIntersection.hh"

class G4Ray;
class G4CylindricalSurface;

class G4SurfaceBoundary
{

public:  // with description


  G4SurfaceBoundary();
  ~G4SurfaceBoundary();
    // Constructor & destructor.

  void Init(const G4CurveVector& bounds0);
    // Initializes with a set of closed curves, each of which is an
    // (inner or outer) boundary. No responsibility to delete the curves
    // is taken.

  inline const G4CurveVector& GetBounds() const;
    // Returns closed curve boundaries.

  inline const G4BoundingBox3D& BBox() const;
    // Returns the bounding-box.

  G4SurfaceBoundary* Project(const G4Transform3D& tr=
                             G4Transform3D::Identity);
    // Projection onto the xy plane after transformation tr.
    // The returned object is allocated dynamically; it is the caller's
    // responsibility to delete it.
    // In case the projection maps a line into a point, 0 is returned.

  G4int IntersectRay2D(const G4Ray& ray);
    // Intersects a 2D boundary with a ray. The ray is projected onto the
    // xy plane. If no intersection 0 is returned, otherwise 1 is returned.
    // and the intersection set to intersection0.
    // The intersection point is: ray.start+ray.dir*intersection0.

  G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);
    // Tangent vector to a curve at the point with parameter u.
    // Returns true if exists. The vector is stored in v.


public:  // without description

  void SplitWithPlane(const G4Point3D& p0, 
		      const G4Vector3D& n,
		      G4SurfaceBoundary*& new1, 
		      G4SurfaceBoundary*& new2 );
    // Splits a boundary with a plane containing p0 with normal n.
    // Pointers to the resulting boundaries are put into new1 and new2.
    // It is the caller's responsibility to delete them.
    // To be implemented yet.

  void SplitWithCylinder(const G4CylindricalSurface& c,
			 G4SurfaceBoundary*& new1, 
			 G4SurfaceBoundary*& new2 );
    // Splits a boundary with a cylindrical surface.
    // Pointers to the resulting boundaries are put into new1 and new2.
    // It is the caller's responsibility to delete them.
    // To be implemented yet.

  inline G4int GetNumberOfPoints() const;
  inline const G4Point3D& GetPoint(G4int Count) const;
    // Functions probably not used and should be removed in the future

  // void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);

public:

  G4Point3DVector points;

private:

  // copy disabled
  G4SurfaceBoundary(const G4SurfaceBoundary&);
  G4SurfaceBoundary& operator=(const G4SurfaceBoundary&);

private:

  G4CurveVector bounds;
  G4BoundingBox3D bBox;

  // to speed up the tangent computation
  G4CurveRayIntersection lastIntersection;

};

#include "G4SurfaceBoundary.icc"

#endif
