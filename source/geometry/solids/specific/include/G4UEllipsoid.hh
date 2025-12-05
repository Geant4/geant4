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
// G4UEllipsoid
//
// Class description:
//
// Wrapper class for G4Ellipsoid to make use of VecGeom Ellipsoid.

// Author: Gabriele Cosmo (CERN), 13.09.2019
// --------------------------------------------------------------------
#ifndef G4UELLIPSOID_HH
#define G4UELLIPSOID_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedEllipsoid.h>

#include "G4Polyhedron.hh"

/**
 * @brief G4UEllipsoid is a wrapper class for G4Ellipsoid to make use
 * of VecGeom Ellipsoid.
 */

class G4UEllipsoid : public G4UAdapter<vecgeom::UnplacedEllipsoid>
{
  using Shape_t = vecgeom::UnplacedEllipsoid;
  using Base_t  = G4UAdapter<vecgeom::UnplacedEllipsoid>;

  public:

    /**
     * Constructs an ellipsoid, given its input parameters.
     *  @param[in] name The solid name.
     *  @param[in] pxSemiAxis Semiaxis in X.
     *  @param[in] pySemiAxis Semiaxis in Y.
     *  @param[in] pzSemiAxis Semiaxis in Z.
     *  @param[in] pzBottomCut Optional lower cut plane level in Z.
     *  @param[in] pzTopCut Optional upper cut plane level in Z.
     */
    G4UEllipsoid(const G4String& name,
                       G4double pxSemiAxis,
                       G4double pySemiAxis,
                       G4double pzSemiAxis,
                       G4double pzBottomCut = 0.0,
                       G4double pzTopCut = 0.0);

    /**
     * Default destructor.
     */
    ~G4UEllipsoid() override = default;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Accessors.
     */
    G4double GetDx() const;
    G4double GetDy() const;
    G4double GetDz() const;
    G4double GetSemiAxisMax (G4int i) const;
    G4double GetZBottomCut() const;
    G4double GetZTopCut() const;

    /**
     * Modifiers.
     */
    void SetSemiAxis (G4double x, G4double y, G4double z);
    void SetZCuts (G4double newzBottomCut, G4double newzTopCut);

    /**
     * Returns the type ID, "G4Ellipsoid" of the solid.
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
    G4UEllipsoid( const G4UEllipsoid &source );
    G4UEllipsoid &operator=( const G4UEllipsoid &source );
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UEllipsoid::GetEntityType() const
{
  return "G4Ellipsoid";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
