#ifndef included_G4CurvePoint
#define included_G4CurvePoint

// A class capable of storing both the parametric and the non-parametric
//   representation of a point on a curve.
// The representation is evaluated lazily for efficiency.

#include "G4Curve.hh"
#include "G4Point3D.hh"

class G4CurvePoint 
{

public:

  G4CurvePoint(G4Curve& c0);

  void Init(G4Curve& c0);

  G4Curve& GetCurve() const;

  void Reset();

  void Reset(G4double u0);

  void Reset(const G4Point3D& p0);

  void Reset(G4double u0, const G4Point3D& p0);

  G4double GetPPoint();

  const G4Point3D& GetPoint();


protected:
  
  G4CurvePoint() { }

  G4Curve*  c;
  
  G4Point3D p;
  G4double  u;

  G4int notComputed;
  static const G4int pFlag;
  static const G4int uFlag;
  static const G4int allFlags;

};

#include "G4CurvePoint.icc"

#endif
