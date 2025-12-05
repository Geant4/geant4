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
// G4UCons
//
// Class description:
//
// Wrapper class for G4Cons to make use of VecGeom Cone.

// Author: G.Cosmo (CERN), 30.10.2013
// --------------------------------------------------------------------
#ifndef G4UCONS_HH
#define G4UCONS_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedCone.h>

#include "G4Polyhedron.hh"

/**
 * @brief G4UCons is a wrapper class for G4Cons to make use of VecGeom Cone.
 */

class G4UCons : public G4UAdapter<vecgeom::GenericUnplacedCone>
{
  using Shape_t = vecgeom::GenericUnplacedCone;
  using Base_t = G4UAdapter<vecgeom::GenericUnplacedCone>;

  public:

    /**
     * Constructs a cone with the given name and dimensions.
     *  @param[in] pName The name of the solid.
     *  @param[in] pRmin1 Inside radius at -fDz.
     *  @param[in] pRmax1 Outside radius at -fDz
     *  @param[in] pRmin2 Inside radius at +fDz.
     *  @param[in] pRmax2 Outside radius at +fDz
     *  @param[in] pDZ Half length in Z.
     *  @param[in] pSPhi Starting angle of the segment in radians.
     *  @param[in] pDPhi Delta angle of the segment in radians.
     */
    G4UCons(const G4String& pName,
                  G4double pRmin1, G4double pRmax1,
                  G4double pRmin2, G4double pRmax2,
                  G4double pDz,
                  G4double pSPhi, G4double pDPhi);

    /**
     * Default destructor.
     */
   ~G4UCons() override = default;

    /**
     * Dispatch method for parameterisation replication mechanism and
     * dimension computation.
     */
    void ComputeDimensions( G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep ) override;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Accessors.
     */
    G4double GetInnerRadiusMinusZ() const;
    G4double GetOuterRadiusMinusZ() const;
    G4double GetInnerRadiusPlusZ()  const;
    G4double GetOuterRadiusPlusZ()  const;
    G4double GetZHalfLength()       const;
    G4double GetStartPhiAngle()     const;
    G4double GetDeltaPhiAngle()     const;
    G4double GetSinStartPhi()       const;
    G4double GetCosStartPhi()       const;
    G4double GetSinEndPhi()         const;
    G4double GetCosEndPhi()         const;
  
    /**
     * Modifiers.
     */
    void SetInnerRadiusMinusZ (G4double Rmin1 );
    void SetOuterRadiusMinusZ (G4double Rmax1 );
    void SetInnerRadiusPlusZ  (G4double Rmin2 );
    void SetOuterRadiusPlusZ  (G4double Rmax2 );
    void SetZHalfLength       (G4double newDz );
    void SetStartPhiAngle     (G4double newSPhi, G4bool trig=true);
    void SetDeltaPhiAngle     (G4double newDPhi);

    /**
     * Returns the type ID, "G4Cons" of the solid.
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
                                 G4double& pMin, G4double& pMax) const override;

    /**
     * Returns a generated polyhedron as graphical representations.
     */
    G4Polyhedron* CreatePolyhedron() const override;

    /**
     * Copy constructor and assignment operator.
     */
    G4UCons(const G4UCons& rhs);
    G4UCons& operator=(const G4UCons& rhs); 
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UCons::GetEntityType() const
{
  return "G4Cons";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
