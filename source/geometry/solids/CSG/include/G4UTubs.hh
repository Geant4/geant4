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
// G4UTubs
//
// Class description:
//
// Wrapper class for G4Tubs to make use of VecGeom Tube.

// Author: G.Cosmo (CERN), 30.10.2013
// --------------------------------------------------------------------
#ifndef G4UTUBS_HH
#define G4UTUBS_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedTube.h>

#include "G4Polyhedron.hh"

/**
 * @brief G4UTubs is a wrapper class for G4Tubs to make use of VecGeom Tube.
 */

class G4UTubs : public G4UAdapter<vecgeom::GenericUnplacedTube>
{
  using Shape_t = vecgeom::GenericUnplacedTube;
  using Base_t = G4UAdapter<vecgeom::GenericUnplacedTube>;

  public:

    /**
     * Constructs a tubs with the given name and dimensions.
     * It checks the input parameters, converting angles so 0<sphi+dpshi<=2_PI
     * if pdphi>2PI then reset it to 2PI.
     *  @param[in] pName The name of the solid.
     *  @param[in] pRMin Inner radius.
     *  @param[in] pRMax Outer radius.
     *  @param[in] pDz Half length in Z.
     *  @param[in] pSPhi Starting phi angle in radians.
     *  @param[in] pDPhi Angle of the segment in radians.
     */
    G4UTubs( const G4String& pName,
                   G4double pRMin,
                   G4double pRMax,
                   G4double pDz,
                   G4double pSPhi,
                   G4double pDPhi );

    /**
     * Default destructor.
     */
    ~G4UTubs() override = default;

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
    G4double GetInnerRadius   () const;
    G4double GetOuterRadius   () const;
    G4double GetZHalfLength   () const;
    G4double GetStartPhiAngle () const;
    G4double GetDeltaPhiAngle () const;
    G4double GetSinStartPhi   () const;
    G4double GetCosStartPhi   () const;
    G4double GetSinEndPhi     () const;
    G4double GetCosEndPhi     () const;

    /**
     * Modifiers.
     */
    void SetInnerRadius   (G4double newRMin);
    void SetOuterRadius   (G4double newRMax);
    void SetZHalfLength   (G4double newDz);
    void SetStartPhiAngle (G4double newSPhi, G4bool trig=true);
    void SetDeltaPhiAngle (G4double newDPhi);

    /**
     * Returns the type ID, "G4Tubs" of the solid.
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
    G4UTubs(const G4UTubs& rhs);
    G4UTubs& operator=(const G4UTubs& rhs);
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UTubs::GetEntityType() const
{
  return "G4Tubs";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
