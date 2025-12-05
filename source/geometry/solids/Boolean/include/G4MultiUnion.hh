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
// G4MultiUnion
//
// Class description:
//
// An instance of "G4MultiUnion" constitutes a grouping of several solids.
// The constituent solids are stored with their respective location in a node
// instance. An instance of "G4MultiUnion" is subsequently composed of one
// or several nodes.

// Author: Marek Gayer (CERN), 19.10.2012 - Original implementation from USolids
//         Gabriele Cosmo (CERN) 06.04.2017 - Adapted implementation in Geant4
//                                            for VecGeom migration
// --------------------------------------------------------------------
#ifndef G4MULTIUNION_HH
#define G4MULTIUNION_HH

#include <vector>

#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4Point3D.hh"
#include "G4Vector3D.hh"
#include "G4SurfBits.hh"
#include "G4Voxelizer.hh"

class G4Polyhedron;

/**
 * @brief An instance of G4MultiUnion constitutes a grouping of several solids.
 * The constituent solids are stored with their respective location in a node
 * instance. An instance of G4MultiUnion is subsequently composed of one or
 * several nodes.
 */

class G4MultiUnion : public G4VSolid
{
  friend class G4Voxelizer;

  public:

    /**
     * Empty default constructor.
     */
    G4MultiUnion();


    /**
     * Constructor assigning a name and initialising components.
     */
    G4MultiUnion(const G4String& name);

    /**
     * Default destructor.
     */
    ~G4MultiUnion() override = default;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4MultiUnion(__void__&);

    /**
     * Methods to build the multiple union by adding nodes (by pointer or ref).
     *  @param[in] solid The solid to be added to the structure.
     *  @param[in] trans The 3D transformation relative to the structure.
     */
    void AddNode(G4VSolid& solid, const G4Transform3D& trans);
    void AddNode(G4VSolid* solid, const G4Transform3D& trans);

    /**
     * Copy constructor and assignment operator.
     */
    G4MultiUnion(const G4MultiUnion& rhs);
    G4MultiUnion& operator=(const G4MultiUnion& rhs);

    /**
     * Accessors to retrieve a transformation or a solid, given an index
     * and the total number of solids in the structure.
     */
    inline const G4Transform3D& GetTransformation(G4int index) const;
    inline G4VSolid* GetSolid(G4int index) const;
    inline G4int GetNumberOfSolids()const;

    /**
     * Returns if the given point "aPoint" is inside or not the solid.
     */
    EInside Inside(const G4ThreeVector& aPoint) const override;

    // Safety methods
    G4double DistanceToIn(const G4ThreeVector& aPoint) const override;
    G4double DistanceToOut(const G4ThreeVector& aPoint) const override;
    inline void SetAccurateSafety(G4bool flag);

    /**
     * Returns the distance along the normalised vector "aDirection" to the
     * shape, from the point at offset "aPoint". If there is no intersection,
     * return kInfinity. The first intersection resulting from leaving a
     * surface/volume is discarded. Hence, it is tolerant of points on
     * the surface of the shape.
     */
    G4double DistanceToIn(const G4ThreeVector& aPoint,
                          const G4ThreeVector& aDirection) const override;

    /**
     * Computes distance from a point presumably inside the solid to the solid
     * surface. Ignores first surface along each axis systematically (for points
     * inside or outside. Early returns zero in case the second surface is
     * behind the starting point.
     * The normal vector to the crossed surface is always filled.
     * In the case the considered point is located inside the G4MultiUnion
     * structure, it acts as follows:
     *  - investigation of the candidates for the passed point
     *  - progressive moving of the point towards the surface, along the
     *    provided direction
     *  - processing of the normal.
     *  @param[in] aPoint The reference point in space.
     *  @param[in] aDirection The normalised direction.
     *  @param[in] calcNorm Flag unused.
     *  @param[out] validNorm Unused.
     *  @param[out] aNormalVector The exiting outwards normal vector (undefined
     *              Magnitude). 
     *  @returns The distance value to exit a volume.
     */
    G4double DistanceToOut(const G4ThreeVector& aPoint,
                           const G4ThreeVector& aDirection,
                           const G4bool calcNorm = false,
                           G4bool* validNorm = nullptr,
                           G4ThreeVector* aNormalVector = nullptr) const override;

    /**
     * Methods to compute the distance to enter/exit a volume, given point and
     * direction, in presence of voxels-based optimisation structure or not.
     */
    G4double DistanceToInNoVoxels(const G4ThreeVector& aPoint,
                                  const G4ThreeVector& aDirection) const;
    G4double DistanceToOutVoxels(const G4ThreeVector& aPoint,
                                 const G4ThreeVector& aDirection,
                                       G4ThreeVector* aNormalVector) const;
    G4double DistanceToOutNoVoxels(const G4ThreeVector& aPoint,
                                   const G4ThreeVector& aDirection,
                                         G4ThreeVector* aNormalVector) const;

    /**
     * Returns the outwards pointing unit normal of the shape for the
     * surface closest to the point at offset "aPoint".
     */
    G4ThreeVector SurfaceNormal(const G4ThreeVector& aPoint) const override;

    /**
     * Determines the bounding box for the considered instance of G4MultiUnion.
     *  @param[in] aAxis The axis along which computing the extent.
     *  @param[out] aMin The minimum bounding limit point.
     *  @param[out] aMax The maximum bounding limit point.
     */
    void Extent(EAxis aAxis, G4double& aMin, G4double& aMax) const;

    /**
     * Computes the bounding limits of the solid.
     *  @param[out] aMin The minimum bounding limit point.
     *  @param[out] aMax The maximum bounding limit point.
     */
    void BoundingLimits(G4ThreeVector& aMin, G4ThreeVector& aMax) const override;

    /**
     * Calculates the minimum and maximum extent of a solid, when under the
     * specified transform, and within the specified limits.
     *  @param[in] pAxis The axis along which compute the extent.
     *  @param[in] pVoxelLimit The limiting space dictated by voxels.
     *  @param[in] pTransform The internal transformation applied to the solid.
     *  @param[out] pMin The minimum extent value.
     *  @param[out] pMax The maximum extent value.
     *  @returns True if the solid is intersected by the extent region.
     */
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                           G4double& pMin, G4double& pMax) const override;

    /**
     * Returns an estimate of the structure capacity or surface area.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    /**
     * Returns the number of solids part of the structure.
     */
    G4int GetNumOfConstituents() const override;

    /**
     * Returns false if any of the solids part of the structure is not faceted.
     */
    G4bool IsFaceted() const override;

    /**
     * Returns a new allocated clone of the multi-union structure.
     */
    G4VSolid* Clone() const override ;

    /**
     * Returns the type ID, "G4MultiUnion" of the solid.
     */
    G4GeometryType GetEntityType() const override { return "G4MultiUnion"; }

    /**
     * Finalises and prepares for use, creating the optimisation structure
     * for all solids in the structure. It must be called once before
     * navigation use.
     */
    void Voxelize();

    /**
     * Returns the xoxelised optimisation structure.
     */
    inline G4Voxelizer& GetVoxels() const;

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo(std::ostream& os) const override;

    /**
     * Returns a point (G4ThreeVector) randomly and uniformly generated
     * on the surface of a solid.
     */
    G4ThreeVector GetPointOnSurface() const override;

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo ( G4VGraphicsScene& scene ) const override ;
    G4Polyhedron* CreatePolyhedron () const override ;
    G4Polyhedron* GetPolyhedron () const override;

  private:

    /**
     * Utility methods for safety and distance computation.
     */
    EInside InsideNoVoxels(const G4ThreeVector& aPoint) const;
    EInside InsideWithExclusion(const G4ThreeVector& aPoint,
                                      G4SurfBits* bits = nullptr) const;
    G4int SafetyFromOutsideNumberNode(const G4ThreeVector& aPoint,
                                            G4double& safety) const;
    G4double DistanceToInCandidates(const G4ThreeVector& aPoint,
                                    const G4ThreeVector& aDirection,
                                          std::vector<G4int>& candidates,
                                          G4SurfBits& bits) const;

    /**
     * Conversion utilities.
     */
    inline G4ThreeVector GetLocalPoint(const G4Transform3D& trans,
                                       const G4ThreeVector& gpoint) const;
    inline G4ThreeVector GetLocalVector(const G4Transform3D& trans,
                                       const G4ThreeVector& gvec) const;
    inline G4ThreeVector GetGlobalPoint(const G4Transform3D& trans,
                                       const G4ThreeVector& lpoint) const;
    inline G4ThreeVector GetGlobalVector(const G4Transform3D& trans,
                                       const G4ThreeVector& lvec) const;
    void TransformLimits(G4ThreeVector& min, G4ThreeVector& max,
                         const G4Transform3D& transformation) const;

  private:

    struct G4MultiUnionSurface
    {
      G4ThreeVector point;
      G4VSolid* solid;
    };

    std::vector<G4VSolid*> fSolids;
    std::vector<G4Transform3D> fTransformObjs;
    G4Voxelizer fVoxels;              // Vozelizer for the solid
    G4double fCubicVolume = 0.0;      // Cubic Volume
    G4double fSurfaceArea = 0.0;      // Surface Area
    G4double kRadTolerance;           // Cached radial tolerance
    mutable G4bool fAccurate = false; // Accurate safety (off by default)

    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;
};

#include "G4MultiUnion.icc"

#endif
