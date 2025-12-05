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
// G4UPolycone
//
// Class description:
//
// Wrapper class for G4Polycone to make use of VecGeom Polycone.

// Author: Gabriele Cosmo (CERN), 31.10.2013
// --------------------------------------------------------------------
#ifndef G4UPOLYCONE_HH
#define G4UPOLYCONE_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedPolycone.h>

#include "G4TwoVector.hh"
#include "G4PolyconeSide.hh"
#include "G4PolyconeHistorical.hh"
#include "G4Polyhedron.hh"

/**
 * @brief G4UPolycone is a wrapper class for G4Polycone to make use
 * of VecGeom Polycone.
 */

class G4UPolycone : public G4UAdapter<vecgeom::GenericUnplacedPolycone>
{
  using Shape_t = vecgeom::GenericUnplacedPolycone;
  using Base_t  = G4UAdapter<vecgeom::GenericUnplacedPolycone>;

  public:

    /**
     * Constructs a polycone shape, given its parameters.
     *  @param[in] name The solid name.
     *  @param[in] phiStart Initial Phi starting angle.
     *  @param[in] phiTotal Total Phi angle.
     *  @param[in] numZPlanes Number of Z planes.
     *  @param[in] zPlane Position of Z planes, with Z in increasing order.
     *  @param[in] rInner Tangent distance to inner surface.
     *  @param[in] rOuter Tangent distance to outer surface.
     */
    G4UPolycone(const G4String& name, 
                      G4double phiStart,    // initial phi starting angle
                      G4double phiTotal,    // total phi angle
                      G4int numZPlanes,     // number of z planes
                const G4double zPlane[],    // position of z planes
                const G4double rInner[],    // tangent distance to inner surface
                const G4double rOuter[]  ); // tangent distance to outer surface

    /**
     * Alternative constructor of a polycone shape, given corners coordinates.
     *  @param[in] name The solid name.
     *  @param[in] phiStart Initial Phi starting angle.
     *  @param[in] phiTotal Total Phi angle.
     *  @param[in] numRZ Number of corners in r,Z space.
     *  @param[in] r r coordinates of corners.
     *  @param[in] z Z coordinates of corners.
     */
    G4UPolycone(const G4String& name, 
                      G4double phiStart,    // initial phi starting angle
                      G4double phiTotal,    // total phi angle
                      G4int    numRZ,       // number corners in r,z space
                const G4double r[],         // r coordinate of these corners
                const G4double z[]       ); // z coordinate of these corners

    /**
     * Default destructor.
     */
    ~G4UPolycone() override = default;
  
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
    G4double GetStartPhi()    const;
    G4double GetDeltaPhi()    const;
    G4double GetEndPhi()      const;
    G4double GetSinStartPhi() const;
    G4double GetCosStartPhi() const;
    G4double GetSinEndPhi()   const;
    G4double GetCosEndPhi()   const;
    G4bool IsOpen()           const;
    G4int  GetNumRZCorner()   const;
    G4PolyconeSideRZ GetCorner(G4int index) const;
    G4PolyconeHistorical* GetOriginalParameters() const;

    /**
     * Modifier.
     */
    void SetOriginalParameters(G4PolyconeHistorical* pars);

    /**
     * Clears all parameters and rebuild the shape, for use in divisions.
     */
    G4bool Reset();

    /**
     * Returns the type ID, "G4Polycone" of the solid.
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
    G4UPolycone( const G4UPolycone& source );
    G4UPolycone& operator=( const G4UPolycone& source );

  private:

    /**
     * Generic initializer, called by all constructors.
     */
    void SetOriginalParameters();

  private:

    G4bool fGenericPcon; // true if created through the 2nd generic constructor
    G4PolyconeHistorical fOriginalParameters; // original input parameters

    G4double wrStart;
    G4double wrDelta;
    std::vector<G4TwoVector> rzcorners;
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UPolycone::GetEntityType() const
{
  return "G4Polycone";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
