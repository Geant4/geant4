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
// G4Cons
//
// Class description:
//
// A G4Cons is, in the general case, a Phi segment of a cone, with
// half-length fDz, inner and outer radii specified at -fDz and +fDz.
// The Phi segment is described by a starting fSPhi angle, and the
// +fDPhi delta angle for the shape.
// If the delta angle is >=2*pi, the shape is treated as continuous
// in Phi.
//
//  Member Data:
//
//  fRmin1  inside radius at  -fDz
//  fRmin2  inside radius at  +fDz
//  fRmax1  outside radius at -fDz
//  fRmax2  outside radius at +fDz
//  fDz  half length in z
//
//  fSPhi  starting angle of the segment in radians
//  fDPhi  delta angle of the segment in radians
//
//  fPhiFullCone   Boolean variable used for indicate the Phi Section
//
// Note:
//    Internally fSPhi & fDPhi are adjusted so that fDPhi<=2PI,
//    and fDPhi+fSPhi<=2PI. This enables simpler comparisons to be
//    made with (say) Phi of a point.

// Author: Paul Kent (CERN), 19.3.1994 - Code converted to tolerant geometry
// --------------------------------------------------------------------
#ifndef G4CONS_HH
#define G4CONS_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UCONS 1
#endif

#if defined(G4GEOM_USE_UCONS)
  #define G4UCons G4Cons
  #include "G4UCons.hh"
#else

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4CSGSolid.hh"
#include "G4Polyhedron.hh"

/**
 * @brief G4Cons is, in the general case, a Phi segment of a cone, with
 * half-length fDz, inner and outer radii specified at -fDz and +fDz.
 * The Phi segment is described by a starting fSPhi angle, and the
 * +fDPhi delta angle for the shape.
 * If the delta angle is >=2*pi, the shape is treated as continuous in Phi.
 */

class G4Cons : public G4CSGSolid
{
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
    G4Cons(const G4String& pName,
                 G4double pRmin1, G4double pRmax1,
                 G4double pRmin2, G4double pRmax2,
                 G4double pDz,
                 G4double pSPhi, G4double pDPhi);

    /**
     * Default destructor.
     */
    ~G4Cons() override = default;

    /**
     * Accessors.
     */
    inline G4double GetInnerRadiusMinusZ() const;
    inline G4double GetOuterRadiusMinusZ() const;
    inline G4double GetInnerRadiusPlusZ()  const;
    inline G4double GetOuterRadiusPlusZ()  const;
    inline G4double GetZHalfLength()       const;
    inline G4double GetStartPhiAngle()     const;
    inline G4double GetDeltaPhiAngle()     const;
    inline G4double GetSinStartPhi()       const;
    inline G4double GetCosStartPhi()       const;
    inline G4double GetSinEndPhi()         const;
    inline G4double GetCosEndPhi()         const;

    /**
     * Modifiers.
     */
    inline void SetInnerRadiusMinusZ (G4double Rmin1 );
    inline void SetOuterRadiusMinusZ (G4double Rmax1 );
    inline void SetInnerRadiusPlusZ  (G4double Rmin2 );
    inline void SetOuterRadiusPlusZ  (G4double Rmax2 );
    inline void SetZHalfLength       (G4double newDz );
    inline void SetStartPhiAngle     (G4double newSPhi, G4bool trig=true);
    inline void SetDeltaPhiAngle     (G4double newDPhi);

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
    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep) override;

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
    G4double DistanceToIn (const G4ThreeVector& p,
                           const G4ThreeVector& v) const override;
    G4double DistanceToIn (const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p,
                           const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;

    /**
     * Returns the type ID, "G4Cons" of the solid.
     */
    G4GeometryType GetEntityType() const override;

    /**
     * Returns a random point located and uniformly distributed on the
     * surface of the solid.
     */
    G4ThreeVector GetPointOnSurface() const override;

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
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo( G4VGraphicsScene& scene ) const override;
    G4Polyhedron* CreatePolyhedron() const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4Cons(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4Cons(const G4Cons& rhs) = default;
    G4Cons& operator=(const G4Cons& rhs);

  private:

    /**
     * Resets relevant values to zero.
     */
    inline void Initialize();

    /**
     * Reset relevant flags and angle values.
     */
    inline void CheckSPhiAngle(G4double sPhi);
    inline void CheckDPhiAngle(G4double dPhi);
    inline void CheckPhiAngles(G4double sPhi, G4double dPhi);

    /**
     * Recomputes relevant trigonometric values and cache them.
     */
    inline void InitializeTrigonometry();

    /**
     * Algorithm for SurfaceNormal() following the original specification
     * for points not on the surface.
     */
    G4ThreeVector ApproxSurfaceNormal(const G4ThreeVector& p) const;

  private:

    /** Radial and angular tolerances. */
    G4double kRadTolerance, kAngTolerance;

    /** Radial and angular dimensions. */
    G4double fRmin1, fRmin2, fRmax1, fRmax2, fDz, fSPhi, fDPhi;

    /** Cached trigonometric values. */
    G4double sinCPhi, cosCPhi, cosHDPhi, cosHDPhiOT, cosHDPhiIT,
             sinSPhi, cosSPhi, sinEPhi, cosEPhi;

    /** Flag for identification of section or full cone. */
    G4bool fPhiFullCone = false;

    /** Cached half tolerance values. */
    G4double halfCarTolerance, halfRadTolerance, halfAngTolerance;
};

#include "G4Cons.icc"

#endif

#endif
