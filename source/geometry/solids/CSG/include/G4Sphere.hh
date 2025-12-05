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
// G4Sphere
//
// Class description:
//
// A G4Sphere is, in the general case, a section of a spherical shell,
// between specified phi and theta angles
//
// The phi and theta segments are described by a starting angle,
// and the +ve delta angle for the shape.
// If the delta angle is >=2*pi, or >=pi the shape is treated as
// continuous in phi or theta respectively.
//
// Theta must lie between 0-pi (incl).
//
// Member Data:
//
//   fRmin  inner radius
//   fRmax  outer radius
//
//   fSPhi  starting angle of the segment in radians
//   fDPhi  delta angle of the segment in radians
//
//   fSTheta  starting angle of the segment in radians
//   fDTheta  delta angle of the segment in radians
//
//
// Note:
//      Internally fSPhi & fDPhi are adjusted so that fDPhi<=2PI,
//      and fDPhi+fSPhi<=2PI. This enables simpler comparisons to be
//      made with (say) Phi of a point.

// Author: Paul Kent (CERN), 28.03.1994 - Code converted to tolerant geometry
// --------------------------------------------------------------------
#ifndef G4SPHERE_HH
#define G4SPHERE_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_USPHERE 1
#endif

#if defined(G4GEOM_USE_USPHERE)
  #define G4USphere G4Sphere
  #include "G4USphere.hh"
#else

#include <CLHEP/Units/PhysicalConstants.h>
#include "G4CSGSolid.hh"
#include "G4Polyhedron.hh"

class G4VisExtent;

/**
 * @brief G4Sphere is, in the general case, a section of a spherical shell,
 * between specified phi and theta angles.
 * The phi and theta segments are described by a starting angle and the +ve
 * delta angle for the shape. If the delta angle is >=2*pi, or >=pi the shape
 * is treated as continuous in phi or theta respectively.
 * Theta must lie between [0..pi].
 */

class G4Sphere : public G4CSGSolid
{
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
    G4Sphere(const G4String& pName,
                   G4double pRmin, G4double pRmax,
                   G4double pSPhi, G4double pDPhi,
                   G4double pSTheta, G4double pDTheta);

    /**
     * Default destructor.
     */
    ~G4Sphere() override = default;

    /**
     * Accessors.
     */
    inline G4double GetInnerRadius    () const;
    inline G4double GetOuterRadius    () const;
    inline G4double GetStartPhiAngle  () const;
    inline G4double GetDeltaPhiAngle  () const;
    inline G4double GetStartThetaAngle() const;
    inline G4double GetDeltaThetaAngle() const;
    inline G4double GetSinStartPhi    () const;
    inline G4double GetCosStartPhi    () const;
    inline G4double GetSinEndPhi      () const;
    inline G4double GetCosEndPhi      () const;
    inline G4double GetSinStartTheta  () const;
    inline G4double GetCosStartTheta  () const;
    inline G4double GetSinEndTheta    () const;
    inline G4double GetCosEndTheta    () const;

    /**
     * Modifiers.
     */
    inline void SetInnerRadius    (G4double newRMin);
    inline void SetOuterRadius    (G4double newRmax);
    inline void SetStartPhiAngle  (G4double newSphi, G4bool trig = true);
    inline void SetDeltaPhiAngle  (G4double newDphi);
    inline void SetStartThetaAngle(G4double newSTheta);
    inline void SetDeltaThetaAngle(G4double newDTheta);

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
     * Returns the type ID, "G4Sphere" of the solid.
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
    G4VisExtent GetExtent() const override;
    void DescribeYourselfTo(G4VGraphicsScene& scene) const override;
    G4Polyhedron* CreatePolyhedron() const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4Sphere(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4Sphere(const G4Sphere& rhs) = default;
    G4Sphere& operator=(const G4Sphere& rhs);

  private:

    /**
     * Resets relevant values to zero.
     */
    inline void Initialize();

    /**
     * Reset relevant flags and angle values.
     */
    inline void CheckThetaAngles(G4double sTheta, G4double dTheta);
    inline void CheckSPhiAngle(G4double sPhi);
    inline void CheckDPhiAngle(G4double dPhi);
    inline void CheckPhiAngles(G4double sPhi, G4double dPhi);

    /**
     * Recompute relevant trigonometric values and cache them.
     */
    inline void InitializePhiTrigonometry();
    inline void InitializeThetaTrigonometry();

    /**
     * Algorithm for SurfaceNormal() following the original specification
     * for points not on the surface.
     */
    G4ThreeVector ApproxSurfaceNormal(const G4ThreeVector& p) const;

  private:

    /** Radial and angular tolerances. */
    G4double fRminTolerance, fRmaxTolerance, kAngTolerance,
             kRadTolerance, fEpsilon = 2.e-11;

    /** Radial and angular dimensions. */
    G4double fRmin, fRmax, fSPhi, fDPhi, fSTheta, fDTheta;

    /** Cached trigonometric values for Phi angle. */
    G4double sinCPhi, cosCPhi, cosHDPhi, cosHDPhiOT, cosHDPhiIT,
             sinSPhi, cosSPhi, sinEPhi, cosEPhi, hDPhi, cPhi, ePhi;

    /** Cached trigonometric values for Theta angle. */
    G4double sinSTheta, cosSTheta, sinETheta, cosETheta,
             tanSTheta, tanSTheta2, tanETheta, tanETheta2, eTheta;

    /** Flags for identification of section, shell or full sphere. */
    G4bool fFullPhiSphere=false, fFullThetaSphere=false, fFullSphere=true;

    /** Cached half tolerance values. */
    G4double halfCarTolerance, halfAngTolerance;
};

#include "G4Sphere.icc"

#endif

#endif
