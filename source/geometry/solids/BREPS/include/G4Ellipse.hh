// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Ellipse.hh,v 1.6 2000-08-28 15:00:32 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Ellipse
//
// Class description:
// 
// Definition of a generic ellipse curve.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __ELLIPTICCURVE_H
#define __ELLIPTICCURVE_H 

#include "G4CircularCurve.hh"

class G4Ellipse : public G4Conic
{

public:  // with description

  G4Ellipse();
  ~G4Ellipse();
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

  inline G4double GetSemiAxis1() const;
  inline G4double GetSemiAxis2() const;
  inline void Init(const G4Axis2Placement3D& position0,
	           G4double semiAxis10, G4double semiAxis20);
    // Get/Set methods for the geometric data.


public:

  inline G4int IntersectRay2D(const G4Ray& ray);
  // void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);

protected:

  void InitBounded();

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
