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
// G4ClippablePolygon
//
// Class description:
//
// A utility class of a polygon that can be clipped by a voxel.

// Author: David C. Williams (UCSC), 1998
// --------------------------------------------------------------------
#ifndef G4CLIPPABLEPOLYGON_HH
#define G4CLIPPABLEPOLYGON_HH

#include <vector>

#include "G4Types.hh"
#include "geomdefs.hh"
#include "G4ThreeVector.hh"

class G4AffineTransform;
class G4VoxelLimits;

/**
 * @brief G4ClippablePolygon in a utility class defining a polygon
 * that can be clipped by a voxel.
 */

class G4ClippablePolygon
{
  using G4ThreeVectorList = std::vector<G4ThreeVector>;

  public:

    /**
     * Default Constructor and Destructor.
     */
    G4ClippablePolygon();
    ~G4ClippablePolygon() = default;
  
    /**
     * Adds a vertex to collection.
     */
    void AddVertexInOrder( const G4ThreeVector& vertex );

    /**
     * Clears the collection of vertices.
     */
    void ClearAllVertices();
  
    /**
     * Accessor and setter for normal vector.
     */
    inline const G4ThreeVector GetNormal() const;
    inline void SetNormal( const G4ThreeVector& newNormal );
  
    /**
     * Clips the polygon along the Cartesian axes, as specified in 'voxelLimit'.
     *  @returns true if the collection of vertices is not empty.
     */
    G4bool Clip( const G4VoxelLimits& voxelLimit );

    /**
     * Clips the polygon while ignoring the indicated axis.
     *  @returns true if the collection of vertices is not empty.
     */
    G4bool PartialClip( const G4VoxelLimits& voxelLimit,
                        const EAxis IgnoreMe );

    /**
     * Clips the polygon along just one axis, as specified in 'voxelLimit'.
     */
    void ClipAlongOneAxis( const G4VoxelLimits& voxelLimit,
                           const EAxis axis );

    /**
     * Computes the polygon extent along the specified 'axis'.
     *  @param[in] axis The Cartesian axis along which computing the extent.
     *  @param[out] min The minimum extent value.
     *  @param[out] max The maximum extent value.
     *  @returns false if invalid polygon (no vertices).
     */
    G4bool GetExtent( const EAxis axis, 
                      G4double& min, G4double& max ) const;

    /**
     * Returns a pointer to the minimum or maximum point along specified 'axis'.
     * Take care! Do not use pointer after destroying parent polygon.
     */
    const G4ThreeVector* GetMinPoint( const EAxis axis ) const;
    const G4ThreeVector* GetMaxPoint( const EAxis axis ) const;

    /**
     * Returns the number of vertices in the polygon.
     */
    inline std::size_t GetNumVertices() const;

    /**
     * Returns true if collection of vertices is empty.
     */
    inline G4bool Empty() const;
  
    /**
     * Decides if the polygon is in "front" of another when viewed along the
     * specified 'axis'. For our purposes here, it is sufficient to use the
     * minimum extent of the polygon along the axis to determine this.
     */
    G4bool InFrontOf(const G4ClippablePolygon& other, EAxis axis) const;

    /**
     * Decides if this polygon is behind another.
     * Remarks in previous method are valid here too.
     */
    G4bool BehindOf(const G4ClippablePolygon& other, EAxis axis) const; 

    /**
     * Gets min/max vertices distance in or out of a plane.
     *  @param[in] pointOnPlane The point on the plane.
     *  @param[in] planeNormal The normal vector to the plane.
     *  @param[out] min The minimum distance from the plane.
     *  @param[out] max The maximum distance from the plane.
     *  @returns false if invalid polygon (no vertices).
     */
    G4bool GetPlanerExtent( const G4ThreeVector& pointOnPlane, 
                            const G4ThreeVector& planeNormal,
                                  G4double& min, G4double& max ) const;

  private:

    /**
     * Clips 'pPolygon' according to 'pVoxelLimits', which must be only
     * limited along one axis, and either the maximum along the axis must be
     * +kInfinity, or the minimum -kInfinity.
     *  @param[in] pPolygon The polygon to clip.
     *  @param[out] outputPolygon The resulting clipped polygon.
     *  @param[in] pVoxelLimit The Cartesian limits.
     */
    void ClipToSimpleLimits( G4ThreeVectorList& pPolygon,
                             G4ThreeVectorList& outputPolygon,
                       const G4VoxelLimits& pVoxelLimit );

  private:

    G4ThreeVectorList vertices;
    G4ThreeVector normal;
    G4double kCarTolerance;
};

#include "G4ClippablePolygon.icc"

#endif
