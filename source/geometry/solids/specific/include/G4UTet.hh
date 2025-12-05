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
// G4UTet
//
// Class description:
//
// Wrapper class for G4Tet to make use of VecGeom Tet.

// Author: Gabriele Cosmo (CERN), 01.11.2013
// --------------------------------------------------------------------
#ifndef G4UTET_HH
#define G4UTET_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedTet.h>

#include "G4Polyhedron.hh"

/**
 * @brief G4UTet is a wrapper class for G4Tet to make use of VecGeom Tet.
 */

class G4UTet : public G4UAdapter<vecgeom::UnplacedTet>
{

  using Shape_t = vecgeom::UnplacedTet;
  using Base_t = G4UAdapter<vecgeom::UnplacedTet>;

  public:

    /**
     * Constructs a tetrahedra, given its parameters.
     *  @param[in] pName The solid name.
     *  @param[in] anchor The anchor point.
     *  @param[in] p2 Point 2.
     *  @param[in] p3 Point 3.
     *  @param[in] p4 Point 4.
     *  @param[in] degeneracyFlag Flag indicating degeneracy of points.
     */
    G4UTet(const G4String& pName,
           const G4ThreeVector& anchor,
           const G4ThreeVector& p2,
           const G4ThreeVector& p3,
           const G4ThreeVector& p4,
                 G4bool* degeneracyFlag = nullptr);

    /**
     * Default destructor.
     */
    ~G4UTet() override = default;

    /**
     * Dispatch method for parameterisation replication mechanism and
     * dimension computation.
     */
    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep) override;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Returns the type ID, "G4Tet" of the solid.
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
    void SetBoundingLimits(const G4ThreeVector& pMin, const G4ThreeVector& pMax);
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
     * Modifier and accessors, for the four vertices of the shape.
     */
    void SetVertices(const G4ThreeVector& anchor,
                     const G4ThreeVector& p1,
                     const G4ThreeVector& p2,
                     const G4ThreeVector& p3,
                     G4bool* degeneracyFlag = nullptr);
    void GetVertices(G4ThreeVector& anchor,
                     G4ThreeVector& p1,
                     G4ThreeVector& p2,
                     G4ThreeVector& p3) const;
    std::vector<G4ThreeVector> GetVertices() const;

    /**
     * Checks if the tetrahedron is degenerate. A tetrahedron is considered
     * as degenerate in case its minimal height is less than the degeneracy
     * tolerance
     *  @returns true if the tetrahedron is degenerate.
     */
    G4bool CheckDegeneracy(const G4ThreeVector& p0,
                           const G4ThreeVector& p1,
                           const G4ThreeVector& p2,
                           const G4ThreeVector& p3) const;

    /**
     * Copy constructor and assignment operator.
     */
    G4UTet(const G4UTet& rhs);
    G4UTet& operator=(const G4UTet& rhs);

  private:

    G4ThreeVector fBmin, fBmax; // bounding box
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UTet::GetEntityType() const
{
  return "G4Tet";
}

inline G4bool G4UTet::IsFaceted() const
{
  return true;
}

#endif  // G4GEOM_USE_USOLIDS

#endif
