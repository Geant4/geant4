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
// G4Tubs
//
// Class description:
//
// A tube or tube segment with curved sides parallel to
// the z-axis. The tube has a specified half-length along
// the z-axis, about which it is centered, and a given
// minimum and maximum radius. A minimum radius of 0
// corresponds to filled tube /cylinder. The tube segment is
// specified by starting and delta angles for phi, with 0
// being the +x axis, PI/2 the +y axis.
// A delta angle of 2PI signifies a complete, unsegmented
// tube/cylinder.
//
// Member Data:
//
//   fRMin  Inner radius
//   fRMax  Outer radius
//   fDz  half length in z
//
//   fSPhi  The starting phi angle in radians,
//          adjusted such that fSPhi+fDPhi<=2PI, fSPhi>-2PI
//
//   fDPhi  Delta angle of the segment.
//
//   fPhiFullTube   Boolean variable used for indicate the Phi Section

// Author: Paul Kent (CERN), 23.01.1994 - First version
// --------------------------------------------------------------------
#ifndef G4TUBS_HH
#define G4TUBS_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UTUBS 1
#endif

#if defined(G4GEOM_USE_UTUBS)
  #define G4UTubs G4Tubs
  #include "G4UTubs.hh"
#else

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4CSGSolid.hh"
#include "G4Polyhedron.hh"

/**
 * @brief G4Tubs is a tube or tube segment with curved sides parallel to
 * the Z-axis. The tube has a specified half-length along the Z-axis, about
 * which it is centered, and a given minimum and maximum radius. A minimum
 * radius of 0 corresponds to filled tube/cylinder. The tube segment is
 * specified by starting and delta angles for phi, with 0 being the +x axis,
 * PI/2 the +y axis. A delta angle of 2PI signifies a complete, unsegmented
 * tube/cylinder.
 */

class G4Tubs : public G4CSGSolid
{
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
    G4Tubs( const G4String& pName,
                  G4double pRMin,
                  G4double pRMax,
                  G4double pDz,
                  G4double pSPhi,
                  G4double pDPhi );

    /**
     * Default destructor.
     */
    ~G4Tubs() override = default;

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
     * Returns the type ID, "G4Tubs" of the solid.
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
    void DescribeYourselfTo (G4VGraphicsScene& scene) const override;
    G4Polyhedron* CreatePolyhedron () const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4Tubs(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4Tubs(const G4Tubs& rhs) = default;
    G4Tubs& operator=(const G4Tubs& rhs);

  protected:

    /**
     * Resets the relevant values to zero.
     */
    inline void Initialize();
 
    /**
     * Methods resetting relevant flags and angle values.
     */
    inline void CheckSPhiAngle(G4double sPhi);
    inline void CheckDPhiAngle(G4double dPhi);
    inline void CheckPhiAngles(G4double sPhi, G4double dPhi);

    /**
     * Recomputes relevant trigonometric values and caches them.
     */
    inline void InitializeTrigonometry();

    /**
     * Computes fast inverse cylindrical (Rxy) radius for points expected to
     * be on a cylindrical surface. Ensures that surface normal vector
     * produced has magnitude with 'normalTolerance' of unit.
     */
    inline G4double FastInverseRxy( const G4ThreeVector& pos, G4double invRad,
                                    G4double normalTolerance ) const;

    /**
     * Algorithm for SurfaceNormal() following the original specification
     * for points not on the surface.
     */
    G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p ) const;

  protected:

    /** Radial and angular tolerances. */
    G4double kRadTolerance, kAngTolerance;

    /** Tolerance of unity for surface normal. */
    static constexpr G4double kNormTolerance = 1.0e-6;

    /** Radial and angular dimensions. */
    G4double fRMin, fRMax, fDz, fSPhi, fDPhi;

    /** Cached trigonometric values. */
    G4double sinCPhi, cosCPhi, cosHDPhi, cosHDPhiOT, cosHDPhiIT,
             sinSPhi, cosSPhi, sinEPhi, cosEPhi;

    /** Flag for identification of section or full tube. */
    G4bool fPhiFullTube;

    /** More cached values - inverse of Rmax, Rmin. */
    G4double fInvRmax, fInvRmin;

    /** Cached half tolerance values. */
    G4double halfCarTolerance, halfRadTolerance, halfAngTolerance;
};

#include "G4Tubs.icc"

#endif

#endif
