// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Curve.hh,v 1.7 2000-11-08 14:22:01 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Curve
//
// Class Description:
//
// Definition of a generic curve.
// To Initialize objects derived from G4Curve:
//   - Construct curve (the constructor takes no parameters)
//   - call Init()
//   - if the curve is bounded, call one of the SetBounds() methods.

// Author: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __CURVE_H
#define __CURVE_H 

#include "geomdefs.hh"
#include "G4Point3D.hh"
#include "G4Vector3D.hh"
#include "G4BoundingBox3D.hh"
#include "G4Transform3D.hh"
#include "G4Ray.hh"

class G4Ray;
class G4CurveRayIntersection;
class G4CurvePoint;
class G4Surface;

class G4Curve
{

 public:  // with description

  G4Curve();
  virtual ~G4Curve();
    // Constructor & destructor.

  G4Curve(const G4Curve& c);
  G4Curve& operator=(const G4Curve& c);
    // Copy contructor and assignment operator.

  inline G4bool operator==(const G4Curve& right) const;
    // Equality operator.

  virtual G4String GetEntityType() const;
    // Returns shape identifier.

  virtual G4Curve* Project(const G4Transform3D& tr=G4Transform3D::Identity)= 0;
    // Projection onto the xy plane after the transformation tr.
    // The returned object is allocated dynamically; it's caller's
    // responsibility to delete it.
    // In case the projection maps two distinct points into one, 0 is returned.
    // NOTE: this should not occur when using projection with
    //       G4SurfaceOfRevolution.

  virtual G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v)= 0;
    // Returns tangent vector v to a curve at the point cp.
    // Returns true if v exists.

  virtual G4int IntersectRay2D(const G4Ray& ray)= 0;
    // Intersects a 2D curve with a ray.
    // The ray is projected onto the xy plane.
    // If no intersection it returns false, otherwise it returns true,
    // the intersection point is ray.start+ray.dir*intersection0.

  inline const G4Point3D& GetStart() const;
  inline const G4Point3D& GetEnd() const;
    // Return start and endpoints in 3D space.

  inline G4double GetPStart() const;
  inline G4double GetPEnd() const;
    // Return start and endpoints in parameter space.

  inline void SetBounds(G4double p1, G4double p2);
  inline void SetBounds(G4double p1, const G4Point3D& p2);
  inline void SetBounds(const G4Point3D& p1, G4double p2);
  inline void SetBounds(const G4Point3D& p1, const G4Point3D& p2);
    // Set start and endpoints, given points as parameter values
    // or 3D points.

  inline G4bool IsBounded() const;
    // Returns if the curve is bounded.

  inline G4bool IsPOn(G4double param) const;
    // Returns if the parameter is on the curve.

  inline void   SetSameSense(G4int sameSense0);
  inline G4int  GetSameSense() const;
    // The sameSense flag can be used to reverse the orientation
    // of the curve (value false).
    // The curves themselves never use the value of this flag;
    // this is just a convenient mean of storing this piece of
    // topological information.

  virtual G4double GetPMax() const = 0;
    // If the parameter space is closed, return the max value
    // if not, return <=0.

  virtual G4Point3D GetPoint(G4double param) const = 0;
    // Return the point in the 3D space, given correspondent parameter.

  virtual G4double GetPPoint(const G4Point3D& p) const = 0;
    // Return parapmeter give a point in space.
    // In case the point is further off the curve than some tolerance
    // the result is undefined.

  const G4BoundingBox3D* BBox() const;
    // Gets the bounding box for the curve
    // If the curve is not bounded, the result is undefined.

  virtual const char* Name() const;
    // Returns type name.

 public:  // without description

  // virtual void Transform(const G4Transform3D& tr);
     // Transformation of the curve.
  // virtual void IntersectRay2D(const G4Ray&, G4CurveRayIntersection&)= 0;

  virtual void SetParentSrfPtr(const G4Surface*);
    // To be moved to a derived class. Is it really needed?

 protected:

  virtual void InitBounded()= 0;
    // This function will be called after the bounds are set.

 protected:

  G4BoundingBox3D bBox;
  G4Point3D start;
  G4Point3D end;
  G4double  pStart;
  G4double  pEnd;
  G4double  pRange;
  G4bool    bounded;
  G4int     sameSense;

 private:
  
  inline void SetStart(const G4Point3D& pt);
  inline void SetStart(G4double p);
  inline void SetEnd(const G4Point3D& p);
  inline void SetEnd(G4double p);
  inline void SetBoundsRest();

};

#include "G4Curve.icc"

#endif
