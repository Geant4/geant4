// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Hyperbola.hh,v 1.6 2000-08-28 15:00:33 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Hyperbola
//
// Class description:
// 
// Definition of a generic hyperbola.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __HYPERBOLICCURVE_H
#define __HYPERBOLICCURVE_H 

#include "G4Conic.hh"

class G4Hyperbola : public G4Conic
{

public:  // with description

  G4Hyperbola();
  ~G4Hyperbola();
    // Constructor & destructor.

  G4Curve* Project(const G4Transform3D& tr=
                   G4Transform3D::Identity);
    // Transforms and projects the curve.

  G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);
    // Returns tangent to curve at a given point, if existing.
    // The tangent is computed from the 3D point representation.

  inline G4double  GetPMax() const;
  inline G4Point3D GetPoint(G4double param) const;
  inline G4double  GetPPoint(const G4Point3D& p) const;
    // Accessors methods.

  inline G4double GetSemiAxis() const;
  inline G4double GetSemiImagAxis() const;
  inline void Init(G4Axis2Placement3D position0,
	           G4double semiAxis0, G4double semiImagAxis0);
    // Get/Set for the geometric data.

public:  // without description

  inline G4int IntersectRay2D(const G4Ray& ray);

  //void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);
  //G4Hyperbola(G4Point3d, G4Point3d, G4Point3d, 
  //            G4Point3d, G4double, G4double );
  //G4Point3d EvaluateByParameterValue(G4double u);
  //G4Point3d GetBoundMax();
  //G4Point3d GetBoundMin();    

protected:

  void InitBounded();

private:

  G4int Inside(const G4Point3D&, const G4Ray&);

  G4Point3D Focus1;
  G4Point3D Focus2;    
  G4Point3D ProjFocus1;
  G4Point3D ProjFocus2; 

  // geometric data
  G4double semiAxis;
  G4double semiImagAxis;

  G4double ratioAxisImagAxis;

  G4Transform3D toUnitHyperbola;

  G4double forTangent;
    // R_1^2/R_2^2
};

#include "G4Hyperbola.icc"

#endif
