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
// G4UOrb
//
// Class description:
//
// Wrapper class for G4Orb to make use of VecGeom Orb.

// Author: G.Cosmo (CERN), 30.10.2013
// --------------------------------------------------------------------
#ifndef G4UORB_HH
#define G4UORB_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedOrb.h>

#include "G4Polyhedron.hh"

/**
 * @brief G4UOrb is a wrapper class for G4Orb to make use of VecGeom Orb.
 */

class G4UOrb : public G4UAdapter<vecgeom::UnplacedOrb>
{
  using Shape_t = vecgeom::UnplacedOrb;
  using Base_t  = G4UAdapter<vecgeom::UnplacedOrb>;

  public:

    /**
     * Constructs a full sphere, given a name and its radius.
     *  @param[in] pName The name of the solid.
     *  @param[in] pRmax Outer radius.
     */
    G4UOrb(const G4String& pName, G4double pRmax);
       
    /**
     * Default destructor.
     */
    ~G4UOrb() override = default;
    
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
     * Accessors and modifiers.
     */
    G4double GetRadius() const;
    void SetRadius(G4double newRmax);
    G4double GetRadialTolerance() const;

    /**
     * Returns the type ID, "G4Orb" of the solid.
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
    G4UOrb(const G4UOrb& rhs);
    G4UOrb& operator=(const G4UOrb& rhs); 
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UOrb::GetEntityType() const
{
  return "G4Orb";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
