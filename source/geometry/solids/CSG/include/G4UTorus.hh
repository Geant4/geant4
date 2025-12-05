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
// G4UTorus
//
// Class description:
//
// Wrapper class for G4Torus to make use of VecGeom Torus.

// Author: Guilherme Lima (FNAL), 19.08.2015
// --------------------------------------------------------------------
#ifndef G4UTORUS_HH
#define G4UTORUS_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedTorus2.h>

#include "G4Polyhedron.hh"

/**
 * @brief G4UTorus is a wrapper class for G4Torus to make use of VecGeom Torus.
 */

class G4UTorus : public G4UAdapter<vecgeom::UnplacedTorus2>
{
  using Shape_t = vecgeom::UnplacedTorus2;
  using Base_t  = G4UAdapter<vecgeom::UnplacedTorus2>;

  public:

    /**
     * Constructs a torus or torus segment with the given name and dimensions.
     *  @param[in] pName The name of the solid.
     *  @param[in] rmin Inner radius.
     *  @param[in] rmax Outer radius.
     *  @param[in] rtor Swept radius of torus.
     *  @param[in] sPhi Starting Phi angle in radians
     *             adjusted such that fSPhi+fDPhi<=2PI, fSPhi>-2PI.
     *  @param[in] dPhi Delta angle of the segment in radians.
     */
    G4UTorus(const G4String& pName,
                   G4double rmin, G4double rmax, G4double rtor,
                   G4double sphi, G4double dphi);

    /**
     * Default destructor.
     */
    ~G4UTorus() override = default;

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
    G4double GetRmin() const;
    G4double GetRmax() const;
    G4double GetRtor() const;
    G4double GetSPhi() const;
    G4double GetDPhi() const;
    G4double GetSinStartPhi() const;
    G4double GetCosStartPhi() const;
    G4double GetSinEndPhi  () const;
    G4double GetCosEndPhi  () const;

    /**
     * Modifiers.
     */
    void SetRmin(G4double arg);
    void SetRmax(G4double arg);
    void SetRtor(G4double arg);
    void SetSPhi(G4double arg);
    void SetDPhi(G4double arg);

    /**
     * Checks and sets all the parameters given in input. Used in constructor.
     */
    void SetAllParameters(G4double arg1, G4double arg2,
                          G4double arg3, G4double arg4, G4double arg5);

    /**
     * Returns the type ID, "G4Torus" of the solid.
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
    G4UTorus(const G4UTorus& rhs);
    G4UTorus& operator=(const G4UTorus& rhs);
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UTorus::GetEntityType() const
{
  return "G4Torus";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
