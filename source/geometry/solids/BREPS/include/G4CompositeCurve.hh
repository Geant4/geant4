// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CompositeCurve.hh,v 1.5 2000-11-08 14:22:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4CompositeCurve
//
// Class description:
// 
// Definition of a generic composition of curves.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4COMPOSITCURVE_H
#define __G4COMPOSITCURVE_H

#include "G4Curve.hh"
#include "G4CurveVector.hh"
#include "G4CurveRayIntersection.hh"
#include "G4Point3DVector.hh"


class G4CompositeCurve : public G4Curve
{

public:  // with description

  G4CompositeCurve();
  virtual ~G4CompositeCurve();
    // Default constructor & destructor

  G4CompositeCurve(const G4Point3DVector& vertices);
    // Constructor creating closed polygons, given the vertices.
    // No call to Init and SetBounds is needed after calling this
    // constructor.

  virtual G4String GetEntityType() const;
    // Returns entity type identifier.

  virtual G4Curve* Project(const G4Transform3D& tr = G4Transform3D::Identity);
    // Project along trasformation tr.
	
  virtual G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);
    // Tangent computed from the 3D point representation.

  virtual G4double  GetPMax() const;
  virtual G4Point3D GetPoint(G4double param) const;
  virtual G4double  GetPPoint(const G4Point3D& p) const;
  const G4CurveVector& GetSegments() const;
    // Acccessors for the geometric data

  void Init(const G4CurveVector& segments0);
    // Initialises geometric data.
    // The object is not responsible for deleting the curves;
    // only a shallow copy of the CurveVector is made.

public:  // without description

  // virtual void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);
  virtual G4int IntersectRay2D(const G4Ray& ray);

protected:  // with description
  
  virtual void InitBounded();
    // Compute bounding box.

private:

  G4CompositeCurve(const G4CompositeCurve&);
  G4CompositeCurve& operator=(const G4CompositeCurve&);
     // Private copy-constructor and assignment operator.

private:
  
  G4CurveVector segments;
  G4CurveRayIntersection lastIntersection;
    // geometric data

};

#include "G4CompositeCurve.icc"

#endif
