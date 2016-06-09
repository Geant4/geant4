//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// ----------------------------------------------------------------------
// Class G4FPlane
//
// Class Description:
//   
// A G4FPlane is a plane created by 3 points or by an origin, an axis and 
// a direction. The plane created is a G4Plane, where his coefficient a, b, 
// c and d are stored. The equation of the plane is:
//                ax + by + cz = d
// 
// This class contain 2 intersection functions :
//       - closest intersection 
//       - intersection by a ray

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, S.Giani, G.Cosmo.
// ----------------------------------------------------------------------
//
// History
// -------
// - SurfaceNormal always returns the direction of NormalX, always containing
//   the correct orientation for all faces (S.Giani).
// - Addition of default argument sense = 1 in the second constructor (S.Giani).
// - The constructor using iVec now properly stores both the internal and
//   external boundaries in the bounds vector (S.Giani).
// - Proper initialization of sameSense in both the constructors (S.Giani). 
// - Addition of third argument (sense) in the second constructor to ensure
//   consistent setting of the normal in all the client code (S.Giani).
// - Proper use of the tolerance in the Intersect function (S.Giani).
// ----------------------------------------------------------------------
#ifndef __PLANESURFACE_H
#define __PLANESURFACE_H

#include "G4Axis2Placement3D.hh"
#include "G4Plane.hh"
#include "G4Surface.hh"


class G4FPlane : public G4Surface
{

public:  // with description

  G4FPlane();
  virtual ~G4FPlane();
    // Default constructor & destructor.

  G4FPlane( const G4Vector3D& direction, 
	    const G4Vector3D& axis     ,	   
	    const G4Point3D&  Pt0,
      G4int sense = 1 );
    // Normal constructor.

  G4FPlane(const G4Point3DVector* pVec, 
	   const G4Point3DVector* iVec= 0,
	   G4int sense = 1);
    // Constructor used by G4BREPSolidBox and G4BREPSolidPolyhedra.

  G4int Intersect(const G4Ray& G4Rayref);
    // Calculates the intersection of the plane and a ray.

  void CalcBBox();
    // Calculates bounding box.

  void Project();
    // Computes the projection of the plane.
  
  inline G4int GetConvex() const;
    // Return plane's convexity, if so.

  inline G4int GetNumberOfPoints() const;
    // Gets the number of the points on the surface boundary.

  inline G4Point3D GetSrfPoint() const;
    // Gets the location point on the surface.

  inline const G4Point3D& GetPoint(G4int Count) const;
    // Gets a surface boundary point.

  void CalcNormal();
    // Computes normal to surface.

  inline G4Vector3D SurfaceNormal(const G4Point3D& Pt) const;
    // Returns normal to surface.

  inline const char* Name() const;
    // Returns the type identifier.

  G4double ClosestDistanceToPoint(const G4Point3D& Pt);
    // Returns the closest distance from point Pt.

  G4double HowNear( const G4Vector3D& x ) const ;
    // Computes the shortest distance from the point x to the G4FPlane.
    // The distance will always be positive.

  inline G4Axis2Placement3D GetPplace() const;
  inline G4Plane GetPplane() const;
    // Accessors to geometrical data.

public:  // without description

  inline G4int MyType() const;
    // Returns the shape type (used in G4BREPSolid).
  
  G4int IsConvex() const; 
    // Returns -1.  (?)

  inline void Deactivate();
    // Deactive, used in G4Surface.

  inline G4Ray* Norm();
    // Returns the normal (used in BREPSolid).

  inline const G4Point3D& GetHitPoint() const;
    // Returns the hit point of the ray on the surface.

protected:

  void InitBounded();

protected:
  
  G4Point3D hitpoint;
    // Hit point of the ray on the surface.

private:

  G4FPlane(const G4FPlane&);
  G4FPlane& operator=(const G4FPlane&);
    // Private copy constructor and assignment operator.

  inline G4int Sign(G4double a) const;

private:

  G4Axis2Placement3D pplace;
  G4Plane            Pl;
  G4Ray*             NormalX;
  G4int              Convex;
  G4SurfaceBoundary* projectedBoundary;

};

#include "G4FPlane.icc"

#endif
