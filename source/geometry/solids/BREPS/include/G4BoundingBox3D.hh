// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BoundingBox3D.hh,v 1.4 2000-11-08 14:22:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BoundingBox3D
//
// Class description:
// 
// Definition of a generic solid's bounding box in the 3D space.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BoundingBox3D_h
#define __G4BoundingBox3D_h 1

#include "G4Ray.hh"
#include "G4Point3D.hh"
#include "G4Vector3D.hh"

class G4BoundingBox3D
{

public:  // with description

  G4BoundingBox3D();
  G4BoundingBox3D(const G4Point3D&);
  G4BoundingBox3D(const G4Point3D&, const G4Point3D&);    
  ~G4BoundingBox3D();
    // Constructors & destructor.

  G4BoundingBox3D(const G4BoundingBox3D& right);
  G4BoundingBox3D& operator=(const G4BoundingBox3D& right);
    // Copy constructor and assignment operator.

  void Init(const G4Point3D&);
  void Init(const G4Point3D&, const G4Point3D&);
  void Extend(const G4Point3D&);
    // To create/extend the bounding box

  G4Point3D GetBoxMin() const;
  G4Point3D GetBoxMax() const;
  G4double GetDistance() const;
  void SetDistance(G4double distance0);
    // Accessors.

  G4int Inside(const G4Point3D&) const;
    // Returns 1 if the point is inside and on the bbox.
    // Returns 0 if the point is outside the bbox.

public:

  G4int GetTestResult() const;
  G4int Test(const G4Ray&);

  static const G4BoundingBox3D space;

private:

  G4int BoxIntersect(const G4Point3D&, 
		     const G4Point3D&, 
		     const G4Vector3D&) const;

  G4double DistanceToIn(const G4Point3D&,
			const G4Vector3D&) const;			  
    
private:

  G4Point3D box_min;
  G4Point3D box_max;
  G4double distance;

  G4int test_result;

  G4Point3D MiddlePoint;
  G4Vector3D GeantBox;    
};

#include "G4BoundingBox3D.icc"

#endif
