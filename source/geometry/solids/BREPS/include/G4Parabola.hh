// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Parabola.hh,v 1.2 1999-01-14 16:01:38 broglia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef __PARABOLICCURVE_H
#define __PARABOLICCURVE_H 

#include "G4Conic.hh"

class G4Parabola : public G4Conic
{
public:
  G4Parabola();
  ~G4Parabola();

  virtual G4Curve* Project(const G4Transform3D& tr = G4Transform3D::Identity);

  virtual G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);

  // virtual void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);
  virtual G4int IntersectRay2D(const G4Ray& ray);

  virtual G4double  GetPMax();
  virtual G4Point3D GetPoint(G4double param);
  virtual G4double  GetPPoint(const G4Point3D& p);

  // STEP
  G4Parabola(STEPentity& Ent);
  G4Parabola(STEPentity& Ent, InstMgr&);      

  //G4Parabola(G4Point3d, G4Point3d, G4double );
  //G4Point3d EvaluateByParameterValue(const G4double u);
  //G4Point3d GetBoundMax();
  //G4Point3d GetBoundMin();   

  // Get/Set for the geometric data
  void Init(const G4Axis2Placement3D& position0, G4double focalDist0);
  double GetFocalDist() const;

protected:

  virtual void InitBounded();

private:

  // geometric data
  double focalDist;

  // for the intersection
  G4Point3D F;
  G4Point3D L0;

};

#include "G4Parabola.icc"

#endif

