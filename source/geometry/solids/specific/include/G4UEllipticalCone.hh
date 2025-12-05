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
// G4UEllipticalCone
//
// Class description:
//
// Wrapper class for G4EllipticalCone to make use of VecGeom EllipticalCone.

// Author: Gabriele Cosmo (CERN), 13.09.2019
// --------------------------------------------------------------------
#ifndef G4UELLIPTICALCONE_HH
#define G4UELLIPTICALCONE_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedEllipticalCone.h>

#include "G4Polyhedron.hh"

/**
 * @brief G4UEllipticalCone is a wrapper class for G4EllipticalCone
 * to make use of VecGeom EllipticalCone.
 */

class G4UEllipticalCone : public G4UAdapter<vecgeom::UnplacedEllipticalCone>
{
  using Shape_t = vecgeom::UnplacedEllipticalCone;
  using Base_t  = G4UAdapter<vecgeom::UnplacedEllipticalCone>;

  public:

    /**
     * Constructs an elliptical cone, with cut in Z.
     *  @param[in] name The solid name.
     *  @param[in] pxSemiAxis Scalar value, defining the scaling along X-axis.
     *  @param[in] pySemiAxis Scalar value, defining the scaling along Y-axis.
     *  @param[in] zMax The Z-coordinate at the apex.
     *  @param[in] pzTopCut Upper cut plane level.
     */
    G4UEllipticalCone(const G4String& name,
                            G4double pxSemiAxis,
                            G4double pySemiAxis,
                            G4double zMax,
                            G4double pzTopCut);

    /**
     * Default destructor.
     */
    ~G4UEllipticalCone() override = default;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Accessors.
     */
    G4double GetSemiAxisMin () const;
    G4double GetSemiAxisMax () const;
    G4double GetSemiAxisX () const;
    G4double GetSemiAxisY () const;
    G4double GetZMax() const;
    G4double GetZTopCut() const;

    /**
     * Modifiers.
     */
    void SetSemiAxis (G4double x, G4double y, G4double z);
    void SetZCut (G4double newzTopCut);

    /**
     * Returns the type ID, "G4EllipticalCone" of the solid.
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
    G4UEllipticalCone( const G4UEllipticalCone& source );
    G4UEllipticalCone& operator=( const G4UEllipticalCone& source );
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UEllipticalCone::GetEntityType() const
{
  return "G4EllipticalCone";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
