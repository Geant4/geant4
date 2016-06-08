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

class G4SurfaceBoundary {

public:

  // Initialize with a set of closed curves, 
  // each of which is an (inner or outer) boundary.
  // no responsibility to delete the curves is taken.
  // shallow copy of G4Curve-s.

  G4SurfaceBoundary();

  void Init(const G4CurveVector& bounds0);

  const G4CurveVector& GetBounds() const { return bounds; }
  
  virtual ~G4SurfaceBoundary();


  // projection onto the xy plane after transformation tr
  // the returned object is allocated dynamically;
  //   it is the caller's responsibility to delete it
  // in case the projection maps a line into a point,
  //   0 is returned

  G4SurfaceBoundary* Project(const G4Transform3D& tr=G4Transform3D::Identity);
  

  // intersect a 2D boundary (probably obtained with Project) with a ray.
  // the ray is projected onto the xy plane.
  // no intersection: return false
  // intersection: return true, and set intersection0
  //   the intersection point is ray.start+ray.dir*intersection0
  //void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);
  G4int IntersectRay2D(const G4Ray& ray);

  // tangent vector to a curve at the point with parameter u
  // true if exists
  // vector comes into v

  G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);


  // split a boundary with a plane containing p0 with normal n.
  // pointers to the resulting boundaries are put into new1 and new2.
  //   it is the caller's responsibility to delete them.

  void SplitWithPlane(const G4Point3D& p0, 
		      const G4Vector3D& n,
		      G4SurfaceBoundary*& new1, 
		      G4SurfaceBoundary*& new2 );



  void SplitWithCylinder(const G4CylindricalSurface& c,
			 G4SurfaceBoundary*& new1, 
			 G4SurfaceBoundary*& new2      );


  const G4BoundingBox3D& BBox() const { return bBox; }

  // the following functions are probably not used
  // and should be removed in the future

  G4Point3DVector points;
  inline int GetNumberOfPoints(){return points.length();}

  inline const G4Point3D& GetPoint(const int Count){return points.ref(Count);}


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

#endif




