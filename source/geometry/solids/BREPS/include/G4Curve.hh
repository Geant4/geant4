// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Curve.hh,v 1.2 1999-01-14 16:01:07 broglia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
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
public:

  // The right way to Initialize objects derived from G4Curve is:
  //   . Construct (the constructor takes no parameters)
  //   . call Init()
  //   . call one of the SetBounds(), if the curve is bounded (most are)

  G4Curve();
  virtual ~G4Curve();

  virtual G4String GetEntityType() const { return "G4Curve"; }

private:

  G4Curve(const G4Curve&);
  G4Curve& operator=(const G4Curve&);

public:

  // transformation of the curve
  //     virtual void Transform(const G4Transform3D& tr);

  // projection onto the xy plane after transformation tr
  // the returned object is allocated dynamically;
  //   it is the caller's responsibility to delete it
  //   in case the projection maps two distinct points into one,
  //   0 is returned
  //   NOTE: this should not occur when using projection
  //     with G4SurfaceOfRevolution.
  //     For other uses this might be too restrictive...
  virtual G4Curve* Project(const G4Transform3D& tr=G4Transform3D::Identity)= 0;
 
  // tangent vector to a curve at the point with parameter u
  // true if exists
  // vector comes into v
  virtual G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v)= 0;


  // intersect a 2D curve (probably obtained with Project) with a ray.
  // the ray is projected onto the xy plane.
  // no intersection: return false
  // intersection: return true, and set intersection0
  //   the intersection point is ray.start+ray.dir*intersection0
  // virtual void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is)= 0;
  virtual G4int IntersectRay2D(const G4Ray& ray)= 0;

  // start and endpoints
  // in 3D space
  const G4Point3D& GetStart() const;
  const G4Point3D& GetEnd() const;


  // in parameter space
  G4double GetPStart() const;
  G4double GetPEnd() const;


  // set start and endpoints
  // four versions, as both points can be given as parameter values
  // or 3D points
  void SetBounds(G4double p1, G4double p2);
  void SetBounds(G4double p1, const G4Point3D& p2);
  void SetBounds(const G4Point3D& p1, G4double p2);
  void SetBounds(const G4Point3D& p1, const G4Point3D& p2);


  // returns if the curve is bounded
  G4bool IsBounded() const;


  // returns if the parameter is on the curve
  G4bool IsPOn(G4double param);


  // the sameSense flag can be used to reverse the orientation
  // of the curve (value false).
  // the curves themselves never use the value of this flag;
  // this is just a convenient means of storing
  // this piece of topological information.
  void   SetSameSense(G4int sameSense0);
  G4int  GetSameSense() const;


  // if the parameter space is closed, return the max value
  // if not, return <=0
  virtual G4double GetPMax()= 0;


  // parameter -> point
  virtual G4Point3D GetPoint(G4double param)= 0;


  // point -> parameter
  // result is undefined
  // if the point is further off the curve than some tolerance
  virtual G4double GetPPoint(const G4Point3D& p)= 0;


  // get the bounding box for the curve
  // this function only works when the curve is bounded!
  // otherwise, the result is undefined.
  const G4BoundingBox3D* BBox() const;

  // To be moved to a derived class
  // really needed?
  virtual void SetParentSrfPtr(const G4Surface* srf){}


  virtual const char* Name(){return "G4Curve";}

  G4bool operator==(const G4Curve& right) const 
  {
    return this == &right;
  }


protected:

  G4BoundingBox3D bBox;

  // This function will be called after the bounds are set:

  virtual void InitBounded()= 0;

private:
  
  void SetStart(const G4Point3D& pt);
  void SetStart(G4double p);
  void SetEnd(const G4Point3D& p);
  void SetEnd(G4double p);
  void SetBoundsRest();

  G4Point3D start;
  G4Point3D end;
  G4double  pStart;
  G4double  pEnd;
  G4double  pRange;
  G4bool    bounded;
  G4int     sameSense;
};


#include "G4Curve.icc"

#endif




