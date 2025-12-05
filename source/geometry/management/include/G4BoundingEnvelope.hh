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
// G4BoundingEnvelope
//
// Class description:
//
// Helper class to facilitate calculation of the extent of a solid
// within the limits defined by a G4VoxelLimits object.
//
// The function CalculateExtent() of a particular solid can create
// a G4BoundingEnvelope object that bounds the solid and then call
// CalculateExtent() of the G4BoundingEnvelope object.
//
// Calculation of extent uses G4Transform3D, thus takes into account
// scaling and reflection, if any.

// Author: Evgueni Tcherniaev (CERN), 25.05.2016 - Initial version
// --------------------------------------------------------------------
#ifndef G4BOUNDINGENVELOPE_HH
#define G4BOUNDINGENVELOPE_HH

#include <vector>
#include "geomdefs.hh"

#include "G4ThreeVector.hh"
#include "G4VoxelLimits.hh"
#include "G4Transform3D.hh"
#include "G4Point3D.hh"
#include "G4Plane3D.hh"

using G4ThreeVectorList = std::vector<G4ThreeVector>;
using G4Polygon3D = std::vector<G4Point3D>;
using G4Segment3D = std::pair<G4Point3D,G4Point3D>;

/**
 * @brief G4BoundingEnvelope is a helper class to facilitate calculation of
 * the extent of a solid within the limits defined by a G4VoxelLimits object.
 * Calculation of the extent takes into account scaling and reflection, if any.
 */

class G4BoundingEnvelope
{
  public:

    /**
     * Constructor from an axis aligned bounding box (AABB).
     *  @param[in] pMin Lower boundary point.
     *  @param[in] pMax Higher boundary point.
     */
    G4BoundingEnvelope(const G4ThreeVector& pMin,
                       const G4ThreeVector& pMax);

    /**
     * Constructor from a sequence of convex polygons, the polygons should have
     * equal numbers of vertices except first and last polygons which may
     * consist of a single vertex.
     *  @param[in] polygons The list of polygons.
     */
    G4BoundingEnvelope(const std::vector<const G4ThreeVectorList*>& polygons);

    /**
     * Constructor from AABB and a sequence of polygons.
     *  @param[in] pMin Lower boundary point.
     *  @param[in] pMax Higher boundary point.
     *  @param[in] polygons The list of polygons.
     */
    G4BoundingEnvelope(const G4ThreeVector& pMin,
                       const G4ThreeVector& pMax,
                       const std::vector<const G4ThreeVectorList*>& polygons);

    /**
     * Default Destructor.
     */
    ~G4BoundingEnvelope() = default;

    /**
     * Analyses the position of the bounding box relative to the voxel.
     * It returns "true" in the case where the value of the extent can be
     * figured out directly from the dimensions of the bounding box, or
     * it is clear that the bounding box and the voxel do not intersect.
     * The reply "false" means that further calculations are needed.
     *  @param[in] pAxis The axis along which compute the extent.
     *  @param[in] pVoxelLimits The limiting space dictated by voxels.
     *  @param[in] pTransform3D Internal transformation applied to the envelope.
     *  @param[out] pMin The minimum extent value.
     *  @param[out] pMax The maximum extent value.
     *  @returns True if the envelope is intersected by the extent region.
     */
    G4bool BoundingBoxVsVoxelLimits(const EAxis pAxis,
                                    const G4VoxelLimits& pVoxelLimits,
                                    const G4Transform3D& pTransform3D,
                                    G4double& pMin, G4double& pMax) const;

    /**
     * Calculates the minimum and maximum extent of the bounding envelope,
     * when under the specified transform, and within the specified limits.
     *  @param[in] pAxis The axis along which compute the extent.
     *  @param[in] pVoxelLimits The limiting space dictated by voxels.
     *  @param[in] pTransform3D Internal transformation applied to the envelope.
     *  @param[out] pMin The minimum extent value.
     *  @param[out] pMax The maximum extent value.
     *  @returns True if the envelope is intersected by the extent region.
     */
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimits,
                           const G4Transform3D& pTransform3D,
                           G4double& pMin, G4double& pMax) const;

  private:

    /**
     * Checks the correctness of the AABB (axis aligned bounding box).
     */
    void CheckBoundingBox();

    /**
     * Checks the correctness of the sequence of convex polygonal bases.
     */
    void CheckBoundingPolygons();

    /**
     * Finds the max scale factor of the transformation.
     */
    G4double FindScaleFactor(const G4Transform3D& pTransform3D) const;

    /**
     * Creates a list of transformed polygons.
     */
    void TransformVertices(const G4Transform3D& pTransform3D,
                                 std::vector<G4Point3D>& pVertices,
                                 std::vector<std::pair<G4int,G4int>>& pBases) const;

    /**
     * Finds the bounding box of a prism.
     */
    void GetPrismAABB(const G4Polygon3D& pBaseA,
                      const G4Polygon3D& pBaseB,
                            G4Segment3D& pAABB) const;

    /**
     * Creates a list of edges of a prism.
     */
    void CreateListOfEdges(const G4Polygon3D& baseA,
                           const G4Polygon3D& baseB,
                                 std::vector<G4Segment3D>& pEdges) const;

    /**
     * Creates a list of planes bounding a prism.
     */
    void CreateListOfPlanes(const G4Polygon3D& baseA,
                            const G4Polygon3D& baseB,
                                  std::vector<G4Plane3D>& pPlanes) const;

    /**
     * Clips a set of edges by G4VoxelLimits.
     */
    G4bool ClipEdgesByVoxel(const std::vector<G4Segment3D>& pEdges,
                            const G4VoxelLimits& pLimits,
                                  G4Segment3D& pExtent) const;

    /**
     * Clips G4VoxelLimits by set of planes bounding a prism.
     */
    void ClipVoxelByPlanes(G4int pBits,
                           const G4VoxelLimits& pLimits,
                           const std::vector<G4Plane3D>& pPlanes,
                           const G4Segment3D& pAABB,
                                 G4Segment3D& pExtent) const;

  private:

    /** Original bounding box. */
    G4ThreeVector fMin, fMax;

    /** Reference to the original sequence of polygonal bases. */
    const std::vector<const G4ThreeVectorList*>* fPolygons = nullptr;
};

#endif
