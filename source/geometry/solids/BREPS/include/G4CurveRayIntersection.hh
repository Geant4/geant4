#ifndef included_G4CurveRayIntersection
#define included_G4CurveRayIntersection

#include "G4CurvePoint.hh"
#include "G4Ray.hh"

class G4CurveRayIntersection: public G4CurvePoint {

// at first, the interface similar to that of G4CurvePoint:

public:

  G4CurveRayIntersection();
  // must be followed by Init!
  // only the distance is set (to infinity)

  G4CurveRayIntersection(G4Curve& c0, const G4Ray& r0);
 
  void Init(G4Curve& c0, const G4Ray& r0);

  const G4Ray& GetRay() const;

  void Reset();

  void ResetPPoint(G4double u0);

  void Reset(const G4Point3D& p0);

  void Reset(G4double u0, const G4Point3D& p0);

  void ResetDistance(G4double d0);

  void Reset(G4double u0, G4double d0);

  void Reset(const G4Point3D& p0, G4double d0);

  void Reset(G4double u0, const G4Point3D& p0, G4double d0);

  G4double GetPPoint();

  const G4Point3D& GetPoint();

  G4double GetDistance();

protected:

  const G4Ray* r;

  G4double d;

  static const G4int dFlag;

// now the additional functionality

public:

  void Update(G4CurveRayIntersection& is);

  void UpdateWithPointOnCurve(G4CurveRayIntersection& is);

};

#include "G4CurveRayIntersection.icc"

#endif
