// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CircularCurve.hh,v 1.5 2000-08-28 15:00:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4CircularCurve
//
// Class description:
// 
// Definition of a generic circle.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __CIRCULARCURVE_H
#define __CIRCULARCURVE_H 

#include "G4Conic.hh"

class G4CircularCurve : public G4Conic
{

public: // with description

  G4CircularCurve();
  virtual ~G4CircularCurve();
    // Constructor & destructor.

  virtual G4Curve* Project(const G4Transform3D& tr=
                           G4Transform3D::Identity);
    // Project along trasformation tr.

  virtual G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);
    // Tangent computed from the 3D point representation as for all conics.
 
  virtual G4double GetPMax() const;
  virtual G4Point3D GetPoint(G4double param) const;
  virtual G4double GetPPoint(const G4Point3D& p) const;
  inline G4double GetRadius() const;
    // Accessors for the geometric data.

  void Init(const G4Axis2Placement3D& position0, G4double radius0);
    // Initialises geometric data.

public: // without description

  // virtual void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);
  virtual G4int IntersectRay2D(const G4Ray& ray);
  
protected:  // with description

  virtual void InitBounded();
    // Compute bounding box.

private:

  G4double radius;
    // Geometric data.

};

#include "G4CircularCurve.icc"

#endif
