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
// G4VSolid
//
// Class description:
//
// Abstract base class for solids, physical shapes that can be tracked through.
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
// clipping of regions for the computation of a solid's extent.
//
// Some visualization/graphics functions are also defined.

// Author: Paul Kent (CERN), 30.06.1995 - Initial version
// --------------------------------------------------------------------
#ifndef G4VSOLID_HH
#define G4VSOLID_HH

#include "G4Types.hh"
#include "G4String.hh"
#include "geomdefs.hh"

class G4AffineTransform;
class G4VoxelLimits;

class G4VPVParameterisation;
class G4VPhysicalVolume;

class G4VGraphicsScene;
class G4Polyhedron;
class G4VisExtent;
class G4DisplacedSolid;

#include "G4ThreeVector.hh"
#include <vector>

using G4ThreeVectorList = std::vector<G4ThreeVector>;
using G4GeometryType = G4String;

/**
 * @brief G4VSolid is an abstract base class for solids, physical shapes that
 * can be tracked through. Each solid has a name, and the constructors and
 * destructors automatically add and subtract them from the G4SolidStore, a
 * singleton 'master' list of available solids.
 */

class G4VSolid
{
  public:

    /**
     * Constructor for G4VSolid. Creates a new shape, with the supplied name.
     * No provision is made for sharing a common name amongst multiple classes.
     *  @param[in] name The solid's name.
     */
    G4VSolid(const G4String& name);

    /**
     * Default Destructor.
     */
    virtual ~G4VSolid();

    /**
     * Copy constructor and assignment operator.
     */
    G4VSolid(const G4VSolid& rhs);
    G4VSolid& operator=(const G4VSolid& rhs);

    /**
     * Equality operator. Returns true only if addresses are the same.
     */
    inline G4bool operator==(const G4VSolid& s) const;

    /**
     * Getter/setter for the shape's name.
     */
    inline G4String GetName() const;
    void SetName(const G4String& name);

    /**
     * Returns the cached geometrical tolerance.
     */
    inline G4double GetTolerance() const;

    /**
     * Computes the bounding limits of the solid.
     *  @param[out] pMin The minimum bounding limit point.
     *  @param[out] pMax The maximum bounding limit point.
     */
    virtual void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    /**
     * Calculates the minimum and maximum extent of the solid, when under the
     * specified transform, and within the specified limits.
     *  @param[in] pAxis The axis along which compute the extent.
     *  @param[in] pVoxelLimit The limiting space dictated by voxels.
     *  @param[in] pTransform The internal transformation applied to the solid.
     *  @param[out] pMin The minimum extent value.
     *  @param[out] pMax The maximum extent value.
     *  @returns True if the solid is intersected by the extent region.
     */
    virtual G4bool CalculateExtent(const EAxis pAxis,
				   const G4VoxelLimits& pVoxelLimit,
				   const G4AffineTransform& pTransform,
				   G4double& pMin, G4double& pMax) const = 0;

    /**
     * Returns the characterisation of a point at offset 'p' respect
     * to the shape.
     *  @param[in] p The point at offset p.
     *  @returns kOutside if the point is outside the shapes boundaries
     *           plus Tolerance/2; kSurface if the point is less than
     *           Tolerance/2 from a surface; kInside otherwise.
     */
    virtual EInside Inside(const G4ThreeVector& p) const = 0;

    /**
     * Returns the outwards pointing unit normal of the shape for the
     * surface closest to the point at offset 'p'.
     *  @param[in] p The point at offset p.
     *  @returns The outwards pointing unit normal.
     */
    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const = 0;

    /**
     * Returns the distance along the normalised vector 'v' to the shape,
     * from the point at offset 'p'. If there is no intersection, returns
     * kInfinity. The first intersection resulting from 'leaving' a
     * surface/volume is discarded. Hence, it is tolerant of points on
     * the surface of the shape.
     *  @param[in] p The point at offset p.
     *  @param[in] v The normalised direction vector.
     *  @returns The distance to enter the shape.
     */
    virtual G4double DistanceToIn(const G4ThreeVector& p,
                                  const G4ThreeVector& v) const = 0;

    /**
     * Calculates the distance to the nearest surface of a shape from an
     * outside point. The distance can be an underestimate.
     *  @param[in] p The point at offset p.
     *  @returns The safety distance to enter the shape.
     */
    virtual G4double DistanceToIn(const G4ThreeVector& p) const = 0;

    /**
     * Returns the distance along the normalised vector 'v' to the shape,
     * from a point at an offset 'p' inside or on the surface of the shape.
     * Intersections with surfaces, when the point is less than Tolerance/2
     * from a surface must be ignored.
     *  @param[in] p The point at offset p.
     *  @param[in] v The normalised direction vector.
     *  @param[in] calcNorm Flag to indicate if to calculate the normal or not.
     *  @param[out] validNorm Flag set to true if the solid lies entirely
     *              behind or on the exiting surface. It is set false if the
     *              solid does not lie entirely behind or on the exiting surface.
     *              'calcNorm' must be true, otherwise it is unused.
     *  @param[out] n The exiting outwards normal vector (undefined Magnitude).
     *              'calcNorm' must be true, otherwise it is unused.
     *  @returns The distance to exit the shape.
     */
    virtual G4double DistanceToOut(const G4ThreeVector& p,
				   const G4ThreeVector& v,
				   const G4bool calcNorm = false,
				   G4bool* validNorm = nullptr,
				   G4ThreeVector* n = nullptr) const = 0;

    /**
     * Calculates the distance to the nearest surface of a shape from an
     * inside point 'p'. The distance can be an underestimate.
     *  @param[in] p The point at offset p.
     *  @returns The safety distance to exit the shape.
     */
    virtual G4double DistanceToOut(const G4ThreeVector& p) const = 0;

    /**
     * Dispatch method for parameterisation replication mechanism and
     * dimension computation. Throws exception if ComputeDimensions() is
     * called from an illegal derived class.
     */
    virtual void ComputeDimensions(G4VPVParameterisation* p,
	                           const G4int n,
                                   const G4VPhysicalVolume* pRep);

    /**
     * Returns an estimation of the solid volume in internal units.
     * This method may be overloaded by derived classes to compute the
     * exact geometrical quantity for solids where this is possible,
     * or anyway to cache the computed value.
     * Note: the computed value is NOT cached.
     */
    virtual G4double GetCubicVolume();

    /**
     * Returns an estimation of the solid surface area in internal units.
     * This method may be overloaded by derived classes to compute the
     * exact geometrical quantity for solids where this is possible,
     * or anyway to cache the computed value.
     * Note: the computed value is NOT cached.
     */
    virtual G4double GetSurfaceArea();

    /**
     * Provides identification of the class of an object
     * (required for persistency).
     */
    virtual G4GeometryType GetEntityType() const = 0;

    /**
     * Returns a random point located on the surface of the solid.
     * Points returned are not necessarily uniformly distributed.
     */
    virtual G4ThreeVector GetPointOnSurface() const;

    /**
     * Returns the number of constituents used for construction of the solid.
     * For non-Boolean solids the return value is one.
     */
    virtual G4int GetNumOfConstituents() const;

    /**
     * Returns true if the solid has only planar faces, false otherwise.
     */
    virtual G4bool IsFaceted() const;

    /**
     * Returns a pointer of a dynamically allocated copy of the solid.
     * Returns a null pointer with warning in case the concrete solid does not
     * implement this method. The caller has responsibility for ownership.
     */
    virtual G4VSolid* Clone() const;

    /**
     * Dumps contents of the solid to a stream.
     */
    virtual std::ostream& StreamInfo(std::ostream& os) const = 0;

    /**
     * Dumps contents of the solid to the standard output.
     */
    inline void DumpInfo() const;

    // Visualization functions

    /**
     * A "double dispatch" function which identifies the solid
     * to the graphics scene for visualization.
     */
    virtual void DescribeYourselfTo (G4VGraphicsScene& scene) const = 0;

    /**
     * Provides extent (bounding box) as possible hint to the graphics view.
     */
    virtual G4VisExtent GetExtent() const;

    /**
     * Creates a Polyhedron used for Visualisation. It is the caller's
     * responsibility to delete it. A null pointer means "not created".
     */
    virtual G4Polyhedron* CreatePolyhedron() const;

    /**
     * Smart access function - creates on request and stores for future
     * access. A null pointer means "not available".
     */
    virtual G4Polyhedron* GetPolyhedron() const;

    /**
     * If the solid is made up from a Boolean operation of two solids,
     * it returns the number 'no' solid. If the solid is not a "Boolean",
     * it returns a null pointer.
     */
    virtual const G4VSolid* GetConstituentSolid(G4int no) const;
    virtual G4VSolid* GetConstituentSolid(G4int no);

    /**
     * If the solid is a "G4DisplacedSolid", it returns a self pointer
     * else it returns a null pointer.
     */
    virtual const G4DisplacedSolid* GetDisplacedSolidPtr() const;
    virtual G4DisplacedSolid* GetDisplacedSolidPtr();

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4VSolid(__void__&);

    /**
     * Calculates the cubic volume only based on the Inside() method.
     * The accuracy is limited by the second argument 'epsilon' or the
     * statistics expressed by 'nStat'.
     *  @param[in] nStat The number of points to generate for the calculation.
     *  @param[in] epsilon The accuracy value.
     */
    G4double EstimateCubicVolume(G4int nStat, G4double epsilon) const;

    /**
     * Calculates the surface area only based on the Inside() method.
     * The accuracy is limited by the second argument 'epsilon' or the
     * statistics expressed by 'nStat'.
     *  @param[in] nStat The number of points to generate for the calculation.
     *  @param[in] epsilon The accuracy value.
     */
    G4double EstimateSurfaceArea(G4int nStat, G4double epsilon) const;

  protected:

    /**
     * Calculates the maximum and minimum extents of the convex polygon
     * 'pPolygon' along the axis 'pAxis', within the limits 'pVoxelLimit'.
     * If the minimum is less than 'pMin', 'pMin' is set to the new minimum.
     * If the maximum is greater than 'pMax', 'pMax' is set to the new maximum.
     * Modifications to 'pPolygon' are made - it is left in an undefined state.
     *  @param[in,out] pPolygon The points defining the convex polygon.
     *  @param[in] pVoxelLimit The limiting space dictated by voxels.
     *  @param[in] pAxis The axis along which compute the extent.
     *  @param[out] pMin The minimum extent value.
     *  @param[out] pMax The maximum extent value.
     */
    void CalculateClippedPolygonExtent(G4ThreeVectorList& pPolygon,
				       const G4VoxelLimits& pVoxelLimit,
				       const EAxis pAxis,
				       G4double& pMin, G4double& pMax) const;

    /**
     * Calculates the maximum and minimum extents of the polygon described
     * by the vertices: pSectionIndex->pSectionIndex+1->
     *                  pSectionIndex+2->pSectionIndex+3->pSectionIndex
     * in the list 'pVertices'.
     * If the minimum is less than 'pMin', 'pMin' is set to the new minimum.
     * If the maximum is greater than 'pMax', 'pMax' is set to the new maximum.
     * No modifications are made to 'pVertices'.
     *  @param[in] pVertices The vertices list defining the convex polygon.
     *  @param[in] pSectionIndex The starting index for vertices.
     *  @param[in] pVoxelLimit The limiting space dictated by voxels.
     *  @param[in] pAxis The axis along which compute the extent.
     *  @param[out] pMin The minimum extent value.
     *  @param[out] pMax The maximum extent value.
     */
    void ClipCrossSection(G4ThreeVectorList* pVertices,
			  const G4int pSectionIndex,
			  const G4VoxelLimits& pVoxelLimit,
			  const EAxis pAxis,
			  G4double& pMin, G4double& pMax) const;

    /**
     * Calculates the maximum and minimum extents of the polygons
     * joining the CrossSections at pSectionIndex->pSectionIndex+3 and
     *                              pSectionIndex+4->pSectionIndex7
     * in the list 'pVertices', within the boundaries of the voxel limits
     * 'pVoxelLimit'.
     * If the minimum is less than 'pMin', 'pMin' is set to the new minimum.
     * If the maximum is greater than 'pMax', 'pMax' is set to the new maximum.
     * No modifications are made to 'pVertices'.
     *  @param[in] pVertices The vertices list defining the convex polygon.
     *  @param[in] pSectionIndex The starting index for vertices.
     *  @param[in] pVoxelLimit The limiting space dictated by voxels.
     *  @param[in] pAxis The axis along which compute the extent.
     *  @param[out] pMin The minimum extent value.
     *  @param[out] pMax The maximum extent value.
     */
    void ClipBetweenSections(G4ThreeVectorList* pVertices,
			     const G4int pSectionIndex,
			     const G4VoxelLimits& pVoxelLimit,
			     const EAxis pAxis,
			     G4double& pMin, G4double& pMax) const;

    /**
     * Clips the specified convex polygon to the given limits, where
     * the polygon is described by the vertices at (0),(1),...,(n),(0) in
     * 'pPolygon'. If the polygon is completely clipped away, the polygon
     * is cleared.
     *  @param[in,out] pPolygon pPolygon The points defining the convex polygon.
     *  @param[in] pVoxelLimit The limiting space dictated by voxels.
     *  @param[in] pAxis The axis along which apply the clipping.
     */
    void ClipPolygon(G4ThreeVectorList& pPolygon,
		     const G4VoxelLimits& pVoxelLimit,
                     const EAxis pAxis) const;

  protected:

    /** Cached geometrical tolerance. */
    G4double kCarTolerance;

  private:

    /**
     * Clips the specified convex polygon to the given limits, storing the
     * result in 'outputPolygon'. The voxel limits must be limited in one
     * *plane* only: this is achieved by having only X or Y or Z limits,
     * and either the minimum or maximum limit set to -+kInfinity respectively.
     *  @param[in,out] pPolygon pPolygon The points defining the convex polygon.
     *  @param[out] outputPolygon The resulting polygon.
     *  @param[in] pVoxelLimit The limiting space dictated by voxels.
     */
    void ClipPolygonToSimpleLimits(G4ThreeVectorList& pPolygon,
				   G4ThreeVectorList& outputPolygon,
			     const G4VoxelLimits& pVoxelLimit) const;

  private:

    /** The shape's name. */
    G4String fshapeName;
};

/// 
/**
 * Streaming operator. Outputs the solid information to the given stream.
 */
std::ostream& operator<<(std::ostream& os, const G4VSolid& e);

#include "G4VSolid.icc"

#endif
