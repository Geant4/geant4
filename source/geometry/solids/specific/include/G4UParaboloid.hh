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
// G4UParaboloid
//
// Class description:
//
// Wrapper class for G4Paraboloid to make use of VecGeom Paraboloid.

// Author: Guilherme Lima (FNAL), 19.08.2015
// --------------------------------------------------------------------
#ifndef G4UPARABOLOID_HH
#define G4UPARABOLOID_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedParaboloid.h>

#include "G4Polyhedron.hh"

/**
 * @brief G4UParaboloid is a wrapper class for G4Paraboloid
 * to make use of VecGeom Paraboloid.
 */

class G4UParaboloid : public G4UAdapter<vecgeom::UnplacedParaboloid>
{
  using Shape_t = vecgeom::UnplacedParaboloid;
  using Base_t  = G4UAdapter<vecgeom::UnplacedParaboloid>;

  public:

    /**
     * Constructs a paraboloid, given its parameters.
     *  @param[in] name The solid name.
     *  @param[in] dz Half length in Z.
     *  @param[in] rlo Radius at -Dz.
     *  @param[in] rhi Radius at +Dz greater than pR1.
     */
    G4UParaboloid(const G4String& name,
                        G4double dz,
                        G4double rlo,
                        G4double rhi);

    /**
     * Default destructor.
     */
    ~G4UParaboloid() override = default;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Accessors.
     */
    G4double GetZHalfLength() const;
    G4double GetRadiusMinusZ() const;
    G4double GetRadiusPlusZ() const;

    /**
     * Modifiers.
     */
    void SetZHalfLength(G4double dz);
    void SetRadiusMinusZ(G4double r1);
    void SetRadiusPlusZ(G4double r2);

    /**
     * Returns the type ID, "G4Paraboloid" of the solid.
     */
    inline G4GeometryType GetEntityType() const override;

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
                           G4double& pmin, G4double& pmax) const override;

    /**
     * Returns a generated polyhedron as graphical representations.
     */
    G4Polyhedron* CreatePolyhedron() const override;

    /**
     * Copy constructor and assignment operator.
     */
    G4UParaboloid( const G4UParaboloid& source );
    G4UParaboloid& operator=( const G4UParaboloid& source );
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UParaboloid::GetEntityType() const
{
  return "G4Paraboloid";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
