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
// G4Paraboloid
//
// Class description:
//
// A G4Paraboloid represents a solid with parabolic profile with possible
// cuts along the Z axis.
//
// Member Data:
//      dz              z half lenght
//      r1              radius at -dz
//      r2              radius at  dz
//      r2  > r1
//
// Equation for the solid:
//      rho^2 <= k1 * z + k2;
//      -dz <= z <= dz
//      r1^2 = k1 * (-dz) + k2
//      r2^2 = k1 * ( dz) + k2

// Author: Lukas Lindroos (CERN), 10.07.2007 - First implementation
// --------------------------------------------------------------------
#ifndef G4PARABOLOID_HH
#define G4PARABOLOID_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UPARABOLOID 1
#endif

#if (defined(G4GEOM_USE_UPARABOLOID) && defined(G4GEOM_USE_SYS_USOLIDS))
  #define G4UParaboloid G4Paraboloid
  #include "G4UParaboloid.hh"
#else

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4VSolid.hh"
#include "G4Polyhedron.hh"

/**
 * @brief G4Paraboloid represents a solid with parabolic profile
 * with possible cuts along the Z axis.
 */

class G4Paraboloid : public G4VSolid
{
  public:

    /**
     * Constructs a paraboloid, given its parameters.
     *  @param[in] pName The solid name.
     *  @param[in] pDz Half length in Z.
     *  @param[in] pR1 Radius at -Dz.
     *  @param[in] pR2 Radius at +Dz greater than pR1.
     */
    G4Paraboloid(const G4String& pName,
                       G4double pDz,
                       G4double pR1,
                       G4double pR2);

    /**
     * Destructor.
     */
    ~G4Paraboloid() override;

    /**
     * Accessors.
     */
    inline G4double GetZHalfLength() const;
    inline G4double GetRadiusMinusZ() const;
    inline G4double GetRadiusPlusZ() const;

    /**
     * Modifiers.
     */
    void SetZHalfLength(G4double dz);
    void SetRadiusMinusZ(G4double R1);
    void SetRadiusPlusZ(G4double R2);

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
     * Concrete implementations of the expected query interfaces for
     * solids, as defined in the base class G4VSolid.
     */
    EInside Inside(const G4ThreeVector& p) const override;
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const override;
    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const override;
    G4double DistanceToIn(const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p,
                           const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;

    /**
     * Returns the type ID, "G4Paraboloid" of the solid.
     */
    G4GeometryType GetEntityType() const override;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo(std::ostream& os) const override;

    /**
     * Returning an estimation of the solid volume (capacity) and
     * surface area, in internal units.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    /**
     * Returns a random point located and uniformly distributed on the
     * surface of the solid.
     */
    G4ThreeVector GetPointOnSurface() const override;

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo(G4VGraphicsScene& scene) const override;
    G4Polyhedron* CreatePolyhedron() const override;
    G4Polyhedron* GetPolyhedron () const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4Paraboloid(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4Paraboloid(const G4Paraboloid& rhs);
    G4Paraboloid& operator=(const G4Paraboloid& rhs);

  private:

    /**
     * Utility method to cache the computation of the solid's surface area.
     */
    G4double CalculateSurfaceArea() const;

  private:

    G4double fSurfaceArea = 0.0;
    G4double fCubicVolume = 0.0;

    // Cached values
    G4double dz = 0.0; // half height
    G4double r1 = 0.0; // radius at -dz
    G4double r2 = 0.0; // radius at  dz
    G4double k1 = 0.0; // k1 = 0.5*(r2*r2 - r1*r1)/dz
    G4double k2 = 0.0; // k2 = 0.5*(r2*r2 + r1*r1)

    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;
};

#include "G4Paraboloid.icc"

#endif  // defined(G4GEOM_USE_UPARABOLOID) && defined(G4GEOM_USE_SYS_USOLIDS)

#endif // G4PARABOLOID_HH
