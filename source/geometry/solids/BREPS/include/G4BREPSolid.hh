// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolid.hh,v 1.8 2000-11-08 14:21:58 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BREPSolid
//
// Class description:
//
// Base class for generic Boundary REPresentation solid.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __SOLID_H
#define __SOLID_H

#include "G4VSolid.hh"	      
#include "G4Surface.hh"
#include "G4Axis2Placement3D.hh"
#include "G4PointRat.hh"	 
#include "G4BoundingBox3D.hh"	 

class G4Ray;

class G4BREPSolid : public G4VSolid
{

public: // with description

  G4BREPSolid(const G4String& name);
  G4BREPSolid(const G4String&, G4Surface**, G4int);
  virtual ~G4BREPSolid();
    // Constructors & destructor

  virtual G4String GetEntityType() const;
    // Returns identifier for solid type entity.
    // A generic BREP solid is considered a "Closed_Shell".

  virtual void Initialize();
    // Calculates the bounding box for solids and surfaces.
    // Converts concave planes to convex.

  G4bool CalculateExtent(const EAxis              pAxis      ,
			 const G4VoxelLimits&     pVoxelLimit,
			 const G4AffineTransform& pTransform ,
			 G4double&                pMin       , 
			 G4double&                pMax        ) const;
    // Calculate the minimum and maximum extent of the solid, when under the
    // specified transform, and within the specified limits. If the solid
    // is not intersected by the region, return false, else return true.

  virtual EInside Inside(register const G4ThreeVector& Pt) const;
    // Determines if the point Pt is inside, outside or on the surface
    // of the solid.

  virtual G4ThreeVector SurfaceNormal(const G4ThreeVector&) const;
    // Calculates the normal of the surface at a point on the surface
    // The sense of the normal depends on the sense of the surface.

  virtual G4double DistanceToIn(const G4ThreeVector&) const;
    // Calculates the shortest distance ("safety") from a point
    // outside the solid to any boundary of this solid.
    // Return 0 if the point is already inside.

  virtual G4double DistanceToIn(register const G4ThreeVector& Pt,
				register const G4ThreeVector& V) const;
    // Calculates the distance from a point Pt outside the solid
    // to the solid's boundary along a specified direction vector V.
    // Note: Intersections with boundaries less than the tolerance must
    //       be ignored if the direction is away from the boundary.

  virtual G4double DistanceToOut(const G4ThreeVector&) const;
    // Calculates the shortest distance ("safety") from a point inside the
    // solid to any boundary of this solid.
    // Return 0 if the point is already outside.	

  virtual G4double DistanceToOut(register const G4ThreeVector& Pt,
				 register const G4ThreeVector& V,
				 const G4bool  calcNorm=false , 
				 G4bool        *validNorm=0   , 
				 G4ThreeVector *n=0             ) const;
    // Calculates the distance from a point inside the solid to the solid`s
    // boundary along a specified direction vector.
    // Return 0 if the point is already outside.
    // Note: If the shortest distance to a boundary is less than the
    //       tolerance, it is ignored. This allows for a point within a
    //       tolerant boundary to leave immediately.

  G4Point3D Scope() const;
    // Utility function to determine the maximum scope of the solid
    // in the coordinates X, Y, Z. Returned as a G4Point3D.

  void DescribeYourselfTo (G4VGraphicsScene& scene) const;
    // Dispatch function which identifies the solid to the graphics scene.

  G4Polyhedron* CreatePolyhedron () const;
  G4NURBS*      CreateNURBS      () const;
    // Create a G4Polyhedron/G4NURBS/...  (It is the caller's responsibility
    // to delete it).  A null pointer means "not created".

  G4int Intersect(register const G4Ray&) const;
    // Gets the roughly calculated closest intersection point for
    // a b_spline and the accurate point for others.

  inline G4Surface* GetSurface(G4int) const;
  inline void Active(G4int) const;  
  inline G4int Active() const;
  inline G4double GetShortestDistance() const;
  inline G4int GetId() const;
  inline void SetId(G4int);
  inline const G4String& GetName() const;
  inline void SetName(const G4String& name);
  inline G4int GetNumberOfFaces() const;
  inline G4int GetNumberOfSolids() const;
  inline const G4Axis2Placement3D* GetPlace() const;
  inline const G4BoundingBox3D*    GetBBox()  const;
    // Accessors methods.

public:

  inline G4double IntersectionDistance() const;
  inline void IntersectionDistance(G4double) const;
    // Gets and sets intersection distance.

  virtual void Reset() const;
    // Resets all distance attributes.

protected:
 
  G4ThreeVectorList* CreateRotatedVertices(const G4AffineTransform&) const;
  G4bool  IsConvex();

  virtual void CalcBBoxes();
    // Calculates the bounding boxes for the surfaces and for the solid.

  void    CheckSurfaceNormals();
  void    RemoveHiddenFaces(register const G4Ray& G4Rayref, G4int) const;   
    // Deactivates the planar faces that are on the "back" side of a solid.
    // B-splines are not handled by this function. Also cases where the ray
    // starting point is Inside the bbox of the solid are ignored as we don't
    // know if the starting point is Inside the actual solid except for
    // axis-oriented box-like solids.

  void    TestSurfaceBBoxes(register const G4Ray&) const;
    // Tests the bounding-box to all surfaces in List.
    // For planar faces the intersection is instead evaluated.

  inline G4int StartInside() const;  
  inline void StartInside(G4int si) const;

  inline void QuickSort( register G4Surface** SrfVec, 
		         register G4int left, register G4int right) const;

protected:

  static G4int        NumberOfSolids;
  static G4Ray        Track;
  static G4double     ShortestDistance;

  G4int               Box, Convex, AxisBox, PlaneSolid;
  G4Axis2Placement3D* place;
  G4BoundingBox3D*    bbox;   
  G4double            intersectionDistance;
  G4int               active;
  G4int               startInside;
  G4int               nb_of_surfaces;
  G4Point3D           intersection_point;
  G4Surface**         SurfaceVec;
  G4double            RealDist;
  G4String            solidname; 
  G4int               Id;
   

private:

  G4BREPSolid(const G4BREPSolid&);
  G4BREPSolid& operator=(const G4BREPSolid&);
    // Private copy constructor and assignment operator.

  G4int IsBox();
  G4int FinalEvaluation(register const G4Ray&, G4int =0) const;

};

#include "G4BREPSolid.icc"

#endif
