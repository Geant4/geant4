// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Ellipse.hh,v 1.4 2000-01-21 13:47:41 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __ELLIPTICCURVE_H
#define __ELLIPTICCURVE_H 

#include "G4CircularCurve.hh"


class G4Ellipse : public G4Conic
{
public:
  G4Ellipse();
  ~G4Ellipse();
 
  virtual G4Curve* Project(const G4Transform3D& tr = G4Transform3D::Identity);

  virtual G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);
  
  //virtual void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);
  virtual G4int IntersectRay2D(const G4Ray& ray);

  virtual G4double  GetPMax();
  virtual G4Point3D GetPoint(G4double param);
  virtual G4double  GetPPoint(const G4Point3D& p);

  // Get/Set for the geometric data
  void Init(const G4Axis2Placement3D& position0,
	    G4double semiAxis10, G4double semiAxis20);
  
  G4double GetSemiAxis1() const;
  G4double GetSemiAxis2() const;


protected:

  virtual void InitBounded();

private:

  // geometric data
  G4double semiAxis1;
  G4double semiAxis2;
  
  G4double ratioAxis2Axis1;
  
  G4Transform3D toUnitCircle;

  G4double forTangent; // -R_1^2/R_2^2
};

#include "G4Ellipse.icc"

#endif
