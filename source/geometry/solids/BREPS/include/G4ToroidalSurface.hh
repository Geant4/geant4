// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ToroidalSurface.hh,v 1.5 2000-11-08 14:22:05 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ToroidalSurface
//
// Class Description:
//   
// Definition of a toroidal surface.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------

#ifndef __G4TOROIDALSURAFCE
#define __G4TOROIDALSURAFCE

#include "G4FPlane.hh"
#include "G4OsloMatrix.hh"

class G4ToroidalSurface : public G4Surface
{

 public:  // with description

  G4ToroidalSurface();
  G4ToroidalSurface(const G4Vector3D&, 
		    const G4Vector3D&,  
		    const G4Vector3D&,  
		    G4double, 
		    G4double);
  virtual ~G4ToroidalSurface();
    // Constructors & destructor.

  inline G4String GetEntityType() const;
    // Returns type identifier of the shape.

  G4int Intersect(const G4Ray&);
    // Determines the intersection of a ray with a torus.	 
    // Returns 1 if the ray intersects the torus.

  void CalcBBox();
    // Computes the bounding-box.

  inline G4Vector3D  	GetDirection() const;
  inline G4Vector3D 	GetAxis()      const;
  inline G4Point3D	GetLocation()  const;
  inline G4double   	GetMinRadius() const;
  inline G4double   	GetMaxRadius() const;
    // Accessors to geometrical data.

  G4double ClosestDistanceToPoint(const G4Point3D& x);
    // Returns the closest distance to the torus surface from a point x.

  virtual G4Vector3D SurfaceNormal(const G4Point3D& Pt) const;
    // Returns the Normal unit vector to the G4ToroidalSurface at a point Pt
    // on (or nearly on) the G4ToroidalSurface.

  inline void MultiplyPointByMatrix(G4Point3D& Base);
  inline void MultiplyVectorByMatrix(G4Vector3D& DCos);
    // Utility methods.

private:

  G4ToroidalSurface(const G4ToroidalSurface&);
  G4ToroidalSurface& operator=(const G4ToroidalSurface&);
    // Private copy constructor and assignment operator.

  inline G4int IsZero(G4double x) const;
  G4int SolveQuartic(G4double c[], G4double s[]);	
  G4int SolveCubic  (G4double c[], G4double s[]);
  G4int SolveQuadric(G4double c[], G4double s[]);
    // Algorithms for solving quadratic, cubic and quartic equations.

private:

  G4Axis2Placement3D Placement;
  G4double           MinRadius;
  G4double	     MaxRadius;
  G4PointMatrix*     TransMatrix;   // transformation matrix  
  G4Point3D          hitpoint;
  const G4double     EQN_EPS;

};

#include "G4ToroidalSurface.icc"

#endif
