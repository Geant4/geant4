// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Surface.hh,v 1.5 2000-08-28 15:00:34 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Surface
//
// Class description:
// 
// Base class for a generic surface.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __surface_h
#define __surface_h 1

#include "geomdefs.hh"
#include "G4CurveVector.hh"
#include "G4PointRat.hh"
#include "G4Ray.hh"
#include "G4BoundingBox3D.hh"
#include "G4STEPEntity.hh"
#include "G4SurfaceBoundary.hh"

class G4Surface : public G4STEPEntity
{

public:  // with description

  G4Surface();
  virtual ~G4Surface();
    // Constructor & destructor.

  G4int operator==( const G4Surface& s );
    // Equality operator.

  virtual G4String GetEntityType() const;
    // Returns type information, needed for STEP output.

  virtual const char* Name() const;
  virtual G4int MyType() const;
    // Returns type information, redundant versions...

  void SetBoundaries(G4CurveVector*);
    // Sets the boundaries of the surface. The curves in the CurveVector
    // must be non-intersecting closed curves.

  virtual G4double HowNear( const G4Vector3D& x ) const; 
    // Returns the distance from the point x to the Surface.
    // The default for a Surface is the distance from the point to the
    // origin. Overridden by derived classes to take into account boundaries.

  virtual G4double ClosestDistanceToPoint(const G4Point3D& Pt);
    // Returns the closest distance to point Pt, as for HowNear() above.
    // This one is used by G4BREPSolid.

  inline G4Vector3D GetOrigin() const;
  inline G4double GetDistance() const;
  inline void SetDistance(G4double Dist);
  inline G4int IsActive() const;
  inline void SetActive(G4int act);
  inline void Deactivate();
  inline void SetSameSense(G4int sameSense0);
  inline G4int GetSameSense() const;
  inline G4BoundingBox3D* GetBBox() const;
    // Get/Set methods for surface's attributes.

  virtual void Reset();
    // Resets basic attributes.

  virtual G4int Intersect(const G4Ray&);
    // Should return the intersection with a ray. Cannot be called from
    // G4Surface base class. Overridden by subclasses.

  virtual G4Vector3D Normal( const G4Vector3D& p ) const;
    // Returns the Normal unit vector to a Surface at the point p on
    // (or nearly on) the Surface. The default is not well defined,
    // so return ( 0, 0, 0 ). Overridden by subclasses.

  virtual void CalcBBox();
    // Calculates the bounds for a bounding box to the surface.
    // The bounding box is used for a preliminary check of intersection.


public:  // without description

  inline static void Project (G4double& Coord, const G4Point3D& Pt, 
			      const G4Plane& Pl);
    // Utility function returning the projection (Coord) of a point Pt
    // on a plane Pl.

  virtual G4double GetUHit() const;
  virtual G4double GetVHit() const;
    // Overriden by BSplineSurface.
    // uhit and vhit are never set.

  virtual G4Point3D Evaluation(const G4Ray& G4Rayref);
  virtual G4int Evaluate(register const G4Ray& Rayref);  
    // For NURBS, there is a two pass intersection algorithm.
    // Sometimes, the result of the cheap one tells us 
    // that execution of the expensive one is not necessary.
    // Evaluation (Evaluate?) is one of them. Better names wanted!

  virtual void Project();
    // Used by BREPSolid. Thus it's probably needed.

  virtual void CalcNormal();
    // Used only in G4FPlane. Should be private to that class?

  virtual G4int IsConvex() const;
    // Only in G4FPlane. BREPSolid::IsConvex uses it.
    // But who uses BREPSolid::IsConvex? Thus: probably not needed.
    // However, knowing if the surface is convex could be used for
    // optimization. 

  virtual G4int GetConvex() const;
    // Only in G4FPlane, but G4BREPSolid uses them.

  virtual G4int GetNumberOfPoints() const;
  virtual const G4Point3D& GetPoint(G4int Count) const;
    // ???

  virtual G4Ray* Norm() const;
  virtual G4Vector3D SurfaceNormal(const G4Point3D& Pt) const = 0;  
    // There is Normal as well -- so what do these do?

/*
  virtual G4double distanceAlongRay( G4int which_way, const G4Ray* ry,
                                     G4Vector3D& p ) const;
    // Returns distance along a Ray (straight line with G4ThreeVec) to leave
    // or enter a Surface. The input variable which_way should be set to +1
    // to indicate leaving a Surface, -1 to indicate entering a Surface.
    // p is the point of intersection of the Ray with the Surface.
    // This is a default function which just gives the distance
    // between the origin of the Ray and the origin of the Surface.
    // Since a generic Surface doesn't have a well-defined Normal, no
    // further checks are Done. It must be overwritten by derived classes.

  virtual G4double G4Surface::distanceAlongHelix( G4int which_way,
                                                  const Helix* hx,
                                                  G4ThreeVec& p ) const;
    // Returns the distance along a Helix to leave or enter a Surface.  
    // The input variable which_way should be set to +1 to indicate
    // leaving a Surface, -1 to indicate entering a Surface.
    // p is the point of intersection of the Helix with the Surface.
    // This is a default function which just gives the distance
    // between the origin of the Helix and the origin of the Surface.
    // Since a generic Surface doesn't have a well-defined Normal, no
    // further checks are Done. It must be overwritten by derived classes.
*/

public:

  G4BoundingBox3D* bbox;
  G4Point3D closest_hit;
  G4Surface* next;

protected:

  virtual void InitBounded();

protected:

  G4SurfaceBoundary surfaceBoundary;
    // The boundaries of the surface.

  G4int Intersected;
    // BSplineSurface and FPlane sets it, no one gets it.

  G4Vector3D origin;
    // Origin of Surface.

  G4int Type;
  G4int AdvancedFace;
  G4int active;
  G4double distance;
  G4double uhit,vhit;
    // Generic attributes...

  G4int sameSense;
    // by L. Broglia

protected:

  const G4double FLT_MAXX;
    // Maybe kInfinity instead?

  const G4double FLT_EPSILO;
    // Maybe kCarTolerance instead?

};

#include "G4Surface.icc"

#endif
