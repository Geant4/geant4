// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VSolid.hh,v 1.2 1999-11-09 18:41:44 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//  Class Description 
//
// class G4VSolid
//
// Abstract base class solids, physical shapes that can be tracked through.
//
// Each solid has a name, and the constructors and destructors automatically
// add and subtract them from the G4SolidStore, a singleton `master' List
// of available solids.
//
// This class defines, but does not implement, functions to compute
// distances to/from the shape. Functions are also defined
// to check whether a point is inside the shape, to return the
// surface normal of the shape at a given point, and to compute
// the extent of the shape. [see descriptions below]
//
// Some protected/private utility functions are implemented for the 
// clipping of regions for the computation of a solid's extent. Note that
// the clipping mechanism is presently inefficient.
// 
// Some visualization/graphics functions are also defined.
//
//
// Member Functions:
//
// G4VSolid(G4String& name)
//   Creates a new shape, with the supplied name
//   No provision is made for sharing a common name amoungst multiple classes.
//   
// G4String GetName() const
//   Returns the current shape's name
// SetName(const G4String& s)
//   Sets the current shape's name
//
// (All remaining functions are pure virtual)
//
// G4bool CalculateExtent(const EAxis pAxis,
//                       const G4VoxelLimit& pVoxelLimit,
//                       const G4AffineTransform& pTransform,
//                       G4double& min, G4double& max)
//
//   Calculate the minimum and maximum extent of the solid, when under the
//   specified transform, and within the specified limits. If the solid does
//   is not intersected by the region, return false, else return true.
//
// EInside Inside(const G4ThreeVector& p)
//   Returns kOutside if the point at offset p is outside the shapes boundaries
//   plus Tolerance/2, kSurface if the point is <=Tolerance/2 from a surface,
//   otherwise kInside.
//
// G4ThreeVector SurfaceNormal(const G4ThreeVector& p)
//   Returns the outwards pointing unit normal of the shape for the
//   surface closest to the point at offset p.
//
// G4double DistanceToIn(const G4ThreeVector& p)
//   Calculate distance to nearest surface of shape from an outside point
//   The distance can be an underestimate.
//
// G4double DistanceToIn(const G4ThreeVector& p, constG4ThreeVector& v)
//   Return distance along the normalised vector v to the shape, from the
//   point at offset p. If there is no intersection, return kInifinity.
//   The first intersection resulting from `leaving' a surface/volume is
//   discarded. Hence, tolerant of points on surface of shape.
//
// G4double DistanceToOut(const G4ThreeVector& p)
//   Calculate distance to nearest surface of shape from an inside point
//   The distance can be an underestimate.
//
// G4double DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
//		          const G4bool calcNorm=false,
//			     G4bool *validNorm=0,G4ThreeVector *n=0;
//   Return distance along the normalised vector v to the shape, from a point
//   at an offset p inside or on the surface of the shape. Intersections with
//   surfaces, when the point is <Tolerance/2 from a surface must be ignored.
//
//   If calcNorm==true:
//       validNorm set true if the solid lies entirely behind or on the
//                 exiting surface.
//                 n set to exiting outwards normal vector(undefined Magnitude)
//       validNorm=false if the solid does not lie entirely behind or on the
//                 exiting surface
//     calcNorm==false:
//       validNorm and n are unused.
//
//   Call as solid.DistanceToOut(p,v) or by specifying all parameters.
//
//
// Type identification 
//                                (required for persistency and STEP interface)
//
//    virtual G4GeometryType  GetEntityType() const = 0;  
//     Provide identification of the class of an object. 
//
//
// Visualization functions:
//
//    virtual void DescribeYourselfTo (G4VGraphicsScene& scene) const = 0;
//      A "double dispatch" function which identifies the solid
//        to the graphics scene.
//    virtual G4VisExtent GetExtent() const = 0;
//      Provides extent (bounding box) as possible hint to graphics view.
//    virtual G4Polyhedron* CreatePolyhedron () const;
//    virtual G4NURBS*      CreateNURBS      () const;
//      Creates a G4Polyhedron/G4NURBS/...  (It is the caller's reponsibility
//        to delete it.)  A null pointer means "not created".
//
//   Class Description - end:
// -----------------------------------------------------------------------------

// Protected functions:
//
// void CalculateClippedPolygonExtent(G4ThreeVectorList& pPolygon,
//				       const G4VoxelLimits& pVoxelLimit,
//				       const EAxis pAxis, 
//				       G4double& pMin, G4double& pMax) const;
//   Calculate the maximum and minimum extents of the convex polygon pPolygon
//   along the axis pAxis, within the limits pVoxelLimit
//
//   If the minimum is <pMin pMin is set to the new minimum
//   If the maximum is >pMax pMax is set to the new maximum
//
//   Modifications to pPolygon are made - it is left in an undefined state
//
//
// void ClipCrossSection(G4ThreeVectorList* pVertices,
//		         const G4int pSectionIndex,
//			 const G4VoxelLimits& pVoxelLimit,
//			 const EAxis pAxis, 
//			 G4double& pMin, G4double& pMax) const;
//
//   Calculate the maximum and minimum extents of the polygon described
//   by the vertices: pSectionIndex->pSectionIndex+1->
//                    pSectionIndex+2->pSectionIndex+3->pSectionIndex
//   in the List pVertices
//
//   If the minimum is <pMin pMin is set to the new minimum
//   If the maximum is >pMax pMax is set to the new maximum
//
//   No modifications are made to pVertices
//
//
// void ClipBetweenSections(G4ThreeVectorList* pVertices,
//			    const G4int pSectionIndex,
//			    const G4VoxelLimits& pVoxelLimit,
//			    const EAxis pAxis, 
//			    G4double& pMin, G4double& pMax) const;
//
//   Calculate the maximum and minimum extents of the polygons
//   joining the CrossSections at pSectionIndex->pSectionIndex+3 and
//                                pSectionIndex+4->pSectionIndex7
//
//   in the List pVertices, within the boundaries of the voxel limits.
//
//   If the minimum is <pMin pMin is set to the new minimum
//   If the maximum is >pMax pMax is set to the new maximum
//
//   No modifications are made to pVertices
//
//
// void ClipPolygon(G4ThreeVectorList& pPolygon,
//	            const G4VoxelLimits& pVoxelLimit) const;
//
//   Clip the specified convex polygon to the given limits, where
//   the polygon is described by the vertices at (0),(1),...,(n),(0) in
//   pPolygon. If the polygon is completely clipped away, the polygon is
//   cleared.
//
//
//
//
// Private functions:
//
// void ClipPolygonToSimpleLimits(G4ThreeVectorList& pPolygon,
//				  G4ThreeVectorList& outputPolygon,
//				  const G4VoxelLimits& pVoxelLimit) const;
//
//   Clip the specified convex polygon to the given limits, storing the
//   result in outputPolygon. The voxel limits must be limited in one
//   *plane* only: This is achieved by having only x or y or z limits,
//   and either the minimum or maximum limit set to -+kInfinity respectively.
//
//
//
//
// Operators:
//
// G4bool operator==(const G4VSolid& s) const
//   Return true only if addresses are the same
//

// Member Data:
//
// G4String fshapeName
//   Name for this solid.
//
// History:
// 17.06.98 J.Apostolakis Added pure virtual function GetEntityType()
// 26.07.96 P.Kent        Added ComputeDimensions for replication mechanism.
// 27.03.96 J.Allison     Methods for visualisation 
// 30.06.95 P.Kent        Initial version, no scoping or visualisation functions

#ifndef G4VSOLID_HH
#define G4VSOLID_HH

#include "globals.hh"
#include "geomdefs.hh"

class G4AffineTransform;
class G4VoxelLimits;

class G4VPVParameterisation;
class G4VPhysicalVolume;

class G4VGraphicsScene;
class G4Polyhedron;
class G4NURBS;
class G4VisExtent;
class G4DisplacedSolid;

#include "G4ThreeVector.hh"
#include <rw/tvordvec.h>
typedef RWTValOrderedVector<G4ThreeVector> G4ThreeVectorList;
typedef G4String   G4GeometryType;

class G4VSolid {
public: // With description 
    G4VSolid(const G4String& name);
    virtual ~G4VSolid();

    G4bool operator==( const G4VSolid& s) const
    {
	return (this==&s) ? true : false;
    }

    G4String GetName() const;
    void SetName(const G4String& name);


    virtual G4bool CalculateExtent(const EAxis pAxis,
				   const G4VoxelLimits& pVoxelLimit,
				   const G4AffineTransform& pTransform,
				   G4double& pMin, G4double& pMax) const = 0;
       
    virtual EInside Inside(const G4ThreeVector& p) const = 0;

    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const = 0;

    virtual G4double DistanceToIn(const G4ThreeVector& p,
                                 const G4ThreeVector& v) const = 0;
    virtual G4double DistanceToIn(const G4ThreeVector& p) const = 0;
    virtual G4double DistanceToOut(const G4ThreeVector& p,
				   const G4ThreeVector& v,
				   const G4bool calcNorm=false,
				   G4bool *validNorm=0,
				   G4ThreeVector *n=0) const = 0;
    virtual G4double DistanceToOut(const G4ThreeVector& p) const = 0;


    virtual void ComputeDimensions(G4VPVParameterisation* p,
	                           const G4int n,
                                   const G4VPhysicalVolume* pRep);
                                   
    virtual G4GeometryType  GetEntityType() const = 0;

    virtual void DescribeYourselfTo (G4VGraphicsScene& scene) const = 0;
    virtual G4VisExtent   GetExtent        () const = 0;
    virtual G4Polyhedron* CreatePolyhedron () const;
    virtual G4NURBS*      CreateNURBS      () const;

    // If Solid is made up from a Boolean operation of two solids,
    //   return the "no" solid. 
    // If the solid is not a "Boolean", return 0
    virtual const G4VSolid* GetConstituentSolid(G4int no) const;
    virtual       G4VSolid* GetConstituentSolid(G4int no);

    // If the Solid is a "G4DisplacedSolid", return a self pointer
    //  else return 0
    virtual const G4DisplacedSolid* GetDisplacedSolidPtr() const; 
    virtual       G4DisplacedSolid* GetDisplacedSolidPtr(); 

protected:

// Calculate the maximum and minimum extents of the convex polygon pPolygon
// along the axis pAxis, within the limits pVoxelLimit
//
// If the minimum is <pMin pMin is set to the new minimum
// If the maximum is >pMax pMax is set to the new maximum
//
// Modifications to pPolygon are made - it is left in an undefined state

    void CalculateClippedPolygonExtent(G4ThreeVectorList& pPolygon,
				       const G4VoxelLimits& pVoxelLimit,
				       const EAxis pAxis, 
				       G4double& pMin, G4double& pMax) const;
// Calculate the maximum and minimum extents of the polygon described
// by the vertices: pSectionIndex->pSectionIndex+1->
//                   pSectionIndex+2->pSectionIndex+3->pSectionIndex
// in the List pVertices
//
// If the minimum is <pMin pMin is set to the new minimum
// If the maximum is >pMax pMax is set to the new maximum
//
// No modifications are made to pVertices
    void ClipCrossSection(G4ThreeVectorList* pVertices,
			  const G4int pSectionIndex,
			  const G4VoxelLimits& pVoxelLimit,
			  const EAxis pAxis, 
			  G4double& pMin, G4double& pMax) const;

// Calculate the maximum and minimum extents of the polygons
// joining the CrossSections at pSectionIndex->pSectionIndex+3 and
//                              pSectionIndex+4->pSectionIndex7
//
// in the List pVertices, within the boundaries of the voxel limits pVoxelLimit
//
// If the minimum is <pMin pMin is set to the new minimum
// If the maximum is >pMax pMax is set to the new maximum
//
// No modifications are made to pVertices

    void ClipBetweenSections(G4ThreeVectorList* pVertices,
			     const G4int pSectionIndex,
			     const G4VoxelLimits& pVoxelLimit,
			     const EAxis pAxis, 
			     G4double& pMin, G4double& pMax) const;

// Clip the specified convex polygon to the given limits, where
// the polygon is described by the vertices at (0),(1),...,(n),(0) in pPolygon.
// 
// If the polygon is completely clipped away, the polygon is cleared.

    void ClipPolygon(G4ThreeVectorList& pPolygon,
		     const G4VoxelLimits& pVoxelLimit) const;

private:

// Clip the specified convex polygon to the given limits, storing the
// result in outputPolygon. The voxel limits must be limited in one
// *plane* only: This is achieved by having only x or y or z limits,
// and either the minimum or maximum limit set to -+kInfinity respectively.

    void ClipPolygonToSimpleLimits(G4ThreeVectorList& pPolygon,
				   G4ThreeVectorList& outputPolygon,
				   const G4VoxelLimits& pVoxelLimit) const;

    G4String fshapeName;	 // Name
};

#include "G4VSolid.icc"

#endif
