#ifndef included_G4CompositeCurve
#define included_G4CompositeCurve

#include "G4Curve.hh"
#include "G4CurveVector.hh"
#include "G4CurveRayIntersection.hh"
#include "G4Point3DVector.hh"


class G4CompositeCurve : public G4Curve
{

public:
  G4CompositeCurve();
  ~G4CompositeCurve();

  // the following is a constructor creating closed polygons,
  // given the vertices.
  // No call to Init and SetBounds is needed after calling this constructor.
  G4CompositeCurve(const G4Point3DVector& vertices);

  virtual G4String GetEntityType() const 
  {
    return "G4CompositeCurve";
  }

  virtual G4Curve* Project(const G4Transform3D& tr = G4Transform3D::Identity);
	
  virtual G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);

  virtual void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);

  virtual G4double  GetPMax();
  virtual G4Point3D GetPoint(G4double param);
  virtual G4double  GetPPoint(const G4Point3D& p);

  // Get/Set for the geometric data
  // the class is not responsible for deleting the curves;
  // only a shallow copy of the CurveVector is made

  void Init(const G4CurveVector& segments0);
  const G4CurveVector& GetSegments() const;


protected:
  
  virtual void InitBounded();
  
private:
  
  // geometric data
  G4CurveVector segments;
  G4CurveRayIntersection lastIntersection;
  
};

#include "G4CompositeCurve.icc"

#endif











