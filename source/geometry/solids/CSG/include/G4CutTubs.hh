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
// G4CutTubs
//
// Class description:
//
// G4CutTubs is a tube with possible cuts in +-Z.
//
// G4CutTubs(pName,pRMin,pRMax,pDZ,pSPhi,pEPhi,pLowNorm,pHighNorm)
//           pName,pRMin,pRMax,pDZ,pSPhi,pEPhi are the same as for G4Tubs,
//           pLowNorm=Outside Normal at -Z
//           pHighNorm=Outside Normal at +Z.

// Author: Tatiana Nikitina (CERN), 31.10.2011
//         Implementation adapted from G4Tubs and TGEo/Ctube implementations.
// --------------------------------------------------------------------
#ifndef G4CUTTUBS_HH
#define G4CUTTUBS_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UCTUBS 1
#endif

#if defined(G4GEOM_USE_UCTUBS)
  #define G4UCutTubs G4CutTubs
  #include "G4UCutTubs.hh"
#else

#include "G4CSGSolid.hh"
#include "G4Polyhedron.hh"

/**
 * @brief G4CutTubs is a tube with possible cuts in +-Z.
 */

class G4CutTubs : public G4CSGSolid
{
  public:

    /**
     * Constructs a tube with the given name and dimensions.
     *  @param[in] pName The name of the solid.
     *  @param[in] pRmin Inner radius.
     *  @param[in] pRmax Outer radius.
     *  @param[in] pDZ Half length in Z.
     *  @param[in] pSPhi Starting angle of the segment in radians.
     *  @param[in] pDPhi Delta angle of the segment in radians.
     *  @param[in] pLowNorm Outside normal vector at -Z.
     *  @param[in] pHighNorm Outside normal vector at +Z.
     */
    G4CutTubs( const G4String& pName,
                     G4double pRMin,
                     G4double pRMax,
                     G4double pDz,
                     G4double pSPhi,
                     G4double pDPhi,
                     G4ThreeVector pLowNorm,
                     G4ThreeVector pHighNorm );

    /**
     * Default destructor.
     */
    ~G4CutTubs() override = default;

    /**
     * Accessors.
     */
    inline G4double GetInnerRadius   () const;
    inline G4double GetOuterRadius   () const;
    inline G4double GetZHalfLength   () const;
    inline G4double GetStartPhiAngle () const;
    inline G4double GetDeltaPhiAngle () const;
    inline G4double GetSinStartPhi   () const;
    inline G4double GetCosStartPhi   () const;
    inline G4double GetSinEndPhi     () const;
    inline G4double GetCosEndPhi     () const;
    inline G4ThreeVector GetLowNorm  () const;
    inline G4ThreeVector GetHighNorm () const;

    /**
     * Modifiers.
     */
    inline void SetInnerRadius   (G4double newRMin);
    inline void SetOuterRadius   (G4double newRMax);
    inline void SetZHalfLength   (G4double newDz);
    inline void SetStartPhiAngle (G4double newSPhi, G4bool trig=true);
    inline void SetDeltaPhiAngle (G4double newDPhi);

    /**
     * Returning an estimation of the solid volume (capacity) and
     * surface area, in internal units.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

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
    EInside Inside( const G4ThreeVector& p ) const override;
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const override;
    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const override;
    G4double DistanceToIn(const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;

    /**
     * Returns the type ID, "G4CutTubs" of the solid.
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
    std::ostream& StreamInfo( std::ostream& os ) const override;

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo ( G4VGraphicsScene& scene ) const override;
    G4Polyhedron* CreatePolyhedron () const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4CutTubs(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4CutTubs(const G4CutTubs& rhs) = default;
    G4CutTubs& operator=(const G4CutTubs& rhs);

  protected:

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
     * Recomputes relevant trigonometric values and caches them.
     */
    inline void InitializeTrigonometry();

    /**
     * Algorithm for SurfaceNormal() following the original specification
     * for points not on the surface.
     */
    G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p ) const;

    /**
     * Checks if the cutted planes are crossing.
     *  @returns True if the solid is ill defined.
     */
    G4bool IsCrossingCutPlanes() const;

    /**
     * Gets the Z value of the point "p" on the cut plane.
     */
    G4double GetCutZ(const G4ThreeVector& p) const;

  private:

    /** Radial and angular tolerances. */
    G4double kRadTolerance, kAngTolerance;

    /** Radial and angular dimensions. */
    G4double fRMin, fRMax, fDz, fSPhi, fDPhi;
    mutable G4double fZMin, fZMax;

    /** Cached trigonometric values. */
    G4double sinCPhi, cosCPhi, cosHDPhi, cosHDPhiOT, cosHDPhiIT,
             sinSPhi, cosSPhi, sinEPhi, cosEPhi;

    /** Flag for identification of section or full tube. */
    G4bool fPhiFullCutTube = false;

    /** Cached half tolerance values. */
    G4double halfCarTolerance, halfRadTolerance, halfAngTolerance;

    /** Normals of Cut at -/+ Dz. */
    G4ThreeVector fLowNorm, fHighNorm;
};

#include "G4CutTubs.icc"

#endif

#endif
