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
// G4UTrd
//
// Class description:
//
// Wrapper class for G4Trd to make use of VecGeom Trd.

// Author: G.Cosmo (CERN), 13.09.2013
// --------------------------------------------------------------------
#ifndef G4UTRD_HH
#define G4UTRD_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedTrd.h>

#include "G4Polyhedron.hh"

/**
 * @brief G4UTrd is a wrapper class for G4Trd to make use of VecGeom Trd.
 */

class G4UTrd : public G4UAdapter<vecgeom::GenericUnplacedTrd> 
{
  using Shape_t = vecgeom::GenericUnplacedTrd;
  using Base_t  = G4UAdapter<vecgeom::GenericUnplacedTrd>;

  public:

    /**
     * Constructs a trapezoid with name, and half lengths.
     *  @param[in] pName The name of the solid.
     *  @param[in] pdx1 Half-length along X at the surface positioned at -dz.
     *  @param[in] pdx2 Half-length along X at the surface positioned at +dz.
     *  @param[in] pdy1 Half-length along Y at the surface positioned at -dz.
     *  @param[in] pdy2 Half-length along Y at the surface positioned at +dz.
     *  @param[in] pdz Half-length along Z axis.
     */
    G4UTrd(const G4String& pName,
                 G4double pdx1, G4double pdx2,
                 G4double pdy1, G4double pdy2,
                 G4double pdz);

    /**
     * Default destructor.
     */
    ~G4UTrd() override = default;

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
     * Accessors.
     */
    G4double GetXHalfLength1() const;
    G4double GetXHalfLength2() const;
    G4double GetYHalfLength1() const;
    G4double GetYHalfLength2() const;
    G4double GetZHalfLength()  const;

    /**
     * Modifiers.
     */
    void SetXHalfLength1(G4double val);
    void SetXHalfLength2(G4double val);
    void SetYHalfLength1(G4double val);
    void SetYHalfLength2(G4double val);
    void SetZHalfLength(G4double val);

    /**
     * Sets all parameters, as for constructor. Checks and sets half-widths.
     */
    void SetAllParameters(G4double pdx1, G4double pdx2,
                          G4double pdy1, G4double pdy2, G4double pdz);

    /**
     * Returns the type ID, "G4Trd" of the solid.
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
    G4UTrd(const G4UTrd& rhs);
    G4UTrd& operator=(const G4UTrd& rhs); 
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UTrd::GetEntityType() const
{
  return "G4Trd";
}

inline G4bool G4UTrd::IsFaceted() const
{
  return true;
}

#endif  // G4GEOM_USE_USOLIDS

#endif
