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
// G4UExtrudedSolid
//
// Class description:
//
// Wrapper class for G4ExtrudedSolid to make use of VecGeom ExtrudedSolid.

// Author: Gabriele Cosmo (CERN), 17.11.2017
// --------------------------------------------------------------------
#ifndef G4UEXTRUDEDSOLID_HH
#define G4UEXTRUDEDSOLID_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedExtruded.h>
#include "G4TwoVector.hh"

#include "G4Polyhedron.hh"

/**
 * @brief G4UExtrudedSolid is a wrapper class for G4ExtrudedSolid
 * to make use of VecGeom ExtrudedSolid.
 */

class G4UExtrudedSolid : public G4UAdapter<vecgeom::UnplacedExtruded>
{
  using Shape_t = vecgeom::UnplacedExtruded;
  using Base_t  = G4UAdapter<vecgeom::UnplacedExtruded>;

  public:

    /**
     * Structure defining a Z section composing the solid.
     */
    struct ZSection
    {
      ZSection() : fZ(0.), fOffset(0.,0.), fScale(1.) {}
      ZSection(G4double z, const G4TwoVector& offset, G4double scale)
        : fZ(z), fOffset(offset), fScale(scale) {}

      G4double    fZ;
      G4TwoVector fOffset;
      G4double    fScale;
    };

    /**
     * General constructor for an extruded polygon, through contour and polyline.
     *  @param[in] pName The solid name.
     *  @param[in] polygon The 2D polygonal contour, i.e. the vertices of the
     *             outlined polygon defined in clock-wise order.
     *  @param[in] zsections The 3D polyline with scale factors, i.e. the
     *             Z-sections defined by Z position in increasing order.
     */
    G4UExtrudedSolid(const G4String&                 pName,
                     const std::vector<G4TwoVector>& polygon,
                     const std::vector<ZSection>&    zsections);

    /**
     * Special constructor for an extruded polygon with 2 Z-sections.
     *  @param[in] pName The solid name.
     *  @param[in] polygon The 2D polygonal contour, i.e. the vertices of the
     *             outlined polygon defined in clock-wise order.
     *  @param[in] halfZ Half length in Z, i.e. the distance from the origin
     *             to the sections.
     *  @param[in] off1 (X, Y) position of the first polygon in -halfZ.
     *  @param[in] scale1 Scale factor at -halfZ.
     *  @param[in] off2 (X, Y) position of the second polygon in +halfZ.
     *  @param[in] scale2 Scale factor at +halfZ.
     */
    G4UExtrudedSolid(const G4String&                 pName,
                     const std::vector<G4TwoVector>& polygon,
                     G4double                        halfZ,
                     const G4TwoVector& off1 = G4TwoVector(0.,0.),
                           G4double scale1 = 1.,
                     const G4TwoVector& off2 = G4TwoVector(0.,0.),
                           G4double scale2 = 1. );

    /**
     * Default Destructor.
     */
    ~G4UExtrudedSolid() override = default;

    /**
     * Accessors.
     */
    G4int GetNofVertices() const;
    G4TwoVector GetVertex(G4int index) const;
    std::vector<G4TwoVector> GetPolygon() const;
    G4int GetNofZSections() const;
    ZSection GetZSection(G4int index) const;
    std::vector<ZSection> GetZSections() const;

    /**
     * Returns the type ID, "G4ExtrudedSolid" of the solid.
     */
    inline G4GeometryType GetEntityType() const override;

    /**
     * Returns true as the solid has only planar faces.
     */
    inline G4bool IsFaceted() const override;

    /**
     * Computes the bounding limits of the solid.
     *  @param[out] pMin The minimum bounding limit point.
     *  @param[out] pMax The maximum bounding limit point.
     */
    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;

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
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const override;  

    /**
     * Returns a generated polyhedron as graphical representations.
     */
    G4Polyhedron* CreatePolyhedron() const override;

    /**
     * Copy constructor and assignment operator.
     */
    G4UExtrudedSolid( const G4UExtrudedSolid& source );
    G4UExtrudedSolid &operator=(const G4UExtrudedSolid& source);
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UExtrudedSolid::GetEntityType() const
{
  return "G4ExtrudedSolid";
}

inline G4bool G4UExtrudedSolid::IsFaceted() const
{
  return true;
}

#endif  // G4GEOM_USE_USOLIDS

#endif
