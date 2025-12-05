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
// G4USphere
//
// Class description:
//
// Wrapper class for G4Sphere to make use of VecGeom Sphere.

// Author: G.Cosmo (CERN), 13.09.2013
// --------------------------------------------------------------------
#ifndef G4USPHERE_HH
#define G4USPHERE_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedSphere.h>

#include "G4Polyhedron.hh"

/**
 * @brief G4USphere is a wrapper class for G4Sphere to make use of
 * VecGeom Sphere.
 */

class G4USphere : public G4UAdapter<vecgeom::UnplacedSphere>
{
  using Shape_t = vecgeom::UnplacedSphere;
  using Base_t  = G4UAdapter<vecgeom::UnplacedSphere>;

  public:

    /**
     * Constructs a sphere or sphere shell section with the given
     * name and dimensions.
     *  @param[in] pName The name of the solid.
     *  @param[in] pRmin Inner radius.
     *  @param[in] pRmax Outer radius.
     *  @param[in] pSPhi Starting Phi angle of the segment in radians.
     *  @param[in] pDPhi Delta Phi angle of the segment in radians.
     *  @param[in] pSTheta Starting Theta angle of the segment in radians.
     *  @param[in] pDTheta Delta Theta angle of the segment in radians.
     */
    G4USphere(const G4String& pName,
                    G4double pRmin, G4double pRmax,
                    G4double pSPhi, G4double pDPhi,
                    G4double pSTheta, G4double pDTheta);
       
    /**
     * Default destructor.
     */
    ~G4USphere() override = default;

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
    G4double GetInnerRadius    () const;
    G4double GetOuterRadius    () const;
    G4double GetStartPhiAngle  () const;
    G4double GetDeltaPhiAngle  () const;
    G4double GetStartThetaAngle() const;
    G4double GetDeltaThetaAngle() const;
    G4double GetSinStartPhi    () const;
    G4double GetCosStartPhi    () const;
    G4double GetSinEndPhi      () const;
    G4double GetCosEndPhi      () const;
    G4double GetSinStartTheta  () const;
    G4double GetCosStartTheta  () const;
    G4double GetSinEndTheta    () const;
    G4double GetCosEndTheta    () const;

    /**
     * Modifiers.
     */
    void SetInnerRadius    (G4double newRMin);
    void SetOuterRadius    (G4double newRmax);
    void SetStartPhiAngle  (G4double newSphi, G4bool trig=true);
    void SetDeltaPhiAngle  (G4double newDphi);
    void SetStartThetaAngle(G4double newSTheta);
    void SetDeltaThetaAngle(G4double newDTheta);

    /**
     * Returns the type ID, "G4Sphere" of the solid.
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
    G4USphere(const G4USphere& rhs);
    G4USphere& operator=(const G4USphere& rhs); 
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4USphere::GetEntityType() const
{
  return "G4Sphere";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
