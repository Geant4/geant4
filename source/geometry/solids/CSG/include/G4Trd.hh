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
// G4Trd
//
// Class description:
//
// A G4Trd is a trapezoid with the x and y dimensions varying along z
// functions:
//
// Member Data:
//
//     fDx1    Half-length along x at the surface positioned at -dz
//     fDx2    Half-length along x at the surface positioned at +dz
//     fDy1    Half-length along y at the surface positioned at -dz
//     fDy2    Half-length along y at the surface positioned at +dz
//     fDz     Half-length along z axis

// Author: Paul Kent (CERN), 12.01.1995 - Code converted to thick geometry
// --------------------------------------------------------------------
#ifndef G4TRD_HH
#define G4TRD_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UTRD 1
#endif

#if defined(G4GEOM_USE_UTRD)
  #define G4UTrd G4Trd
  #include "G4UTrd.hh"
#else

#include "G4CSGSolid.hh"
#include "G4Polyhedron.hh"

/**
 * @brief G4Trd is a trapezoid with the X and Y dimensions varying along Z.
 */

class G4Trd : public G4CSGSolid
{
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
    G4Trd( const G4String& pName,
                 G4double pdx1, G4double pdx2,
                 G4double pdy1, G4double pdy2,
                 G4double pdz );

    /**
     * Default destructor.
     */
   ~G4Trd() override = default;

    /**
     * Accessors.
     */
    inline G4double GetXHalfLength1() const;
    inline G4double GetXHalfLength2() const;
    inline G4double GetYHalfLength1() const;
    inline G4double GetYHalfLength2() const;
    inline G4double GetZHalfLength()  const;

    /**
     * Modifiers.
     */
    inline void SetXHalfLength1(G4double val);
    inline void SetXHalfLength2(G4double val);
    inline void SetYHalfLength1(G4double val);
    inline void SetYHalfLength2(G4double val);
    inline void SetZHalfLength(G4double val);

    /**
     * Sets all parameters, as for constructor. Checks and sets half-widths.
     */
    void SetAllParameters ( G4double pdx1, G4double pdx2,
                            G4double pdy1, G4double pdy2,
                            G4double pdz );

    /**
     * Returning an estimation of the solid volume (capacity) and
     * surface area, in internal units.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    /**
     * Dispatch method for parameterisation replication mechanism and
     * dimension computation.
     */
    void ComputeDimensions( G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep ) override;

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
     * Concrete implementations of the expected query interfaces for
     * solids, as defined in the base class G4VSolid.
     */
    EInside Inside( const G4ThreeVector& p ) const override;
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const override;
    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v ) const override;
    G4double DistanceToIn( const G4ThreeVector& p ) const override;
    G4double DistanceToOut( const G4ThreeVector& p,
                            const G4ThreeVector& v,
                            const G4bool calcNorm = false,
                                  G4bool* validNorm = nullptr,
                                  G4ThreeVector* n = nullptr ) const override;
    G4double DistanceToOut( const G4ThreeVector& p ) const override;

    /**
     * Returns the type ID, "G4Trd" of the solid.
     */
    G4GeometryType GetEntityType() const override;

    /**
     * Returns a random point located and uniformly distributed on the
     * surface of the solid.
     */
    G4ThreeVector GetPointOnSurface() const override;

    /**
     * Returns true as the solid has only planar faces.
     */
    G4bool IsFaceted() const override;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo( std::ostream& os ) const override;

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo (G4VGraphicsScene& scene) const override;
    G4Polyhedron* CreatePolyhedron () const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4Trd(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4Trd(const G4Trd& rhs);
    G4Trd& operator=(const G4Trd& rhs);

  private:

    /**
     * Checks the input parameters.
     */
    void CheckParameters();

    /**
     * Sets the side planes.
     */
    void MakePlanes();

    /**
     * Algorithm for SurfaceNormal() following the original specification
     * for points not on the surface.
     */
    G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p ) const;

  private:

    G4double halfCarTolerance;
    G4double fDx1,fDx2,fDy1,fDy2,fDz,fHx,fHy;
    struct { G4double a,b,c,d; } fPlanes[4];
};

#include "G4Trd.icc"

#endif

#endif
