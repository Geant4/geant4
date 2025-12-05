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
// G4Hype
//
// Class description:
//
// This class implements a tube with hyperbolic profile.
// It describes an hyperbolic volume with curved sides parallel to
// the z-axis. The solid has a specified half-length along the z axis,
// about which it is centered, and a given minimum and maximum radius.
// A minimum radius of 0 signifies a filled Hype (with hyperbolical
// inner surface). To have a filled Hype the user must specify
// inner radius = 0 AND inner stereo angle = 0.
// The inner and outer hyperbolical surfaces can have different
// stereo angles. A stereo angle of 0 gives a cylindrical surface.

// Authors: Ernesto Lamanna & Francesco Safai Tehrani - Created
//          Rome, INFN & University of Rome "La Sapienza", 09.06.1998.
// --------------------------------------------------------------------
#ifndef G4HYPE_HH
#define G4HYPE_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UHYPE 1
#endif

#if (defined(G4GEOM_USE_UHYPE) && defined(G4GEOM_USE_SYS_USOLIDS))
  #define G4UHype G4Hype
  #include "G4UHype.hh"
#else

#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Polyhedron.hh"

class G4SolidExtentList;
class G4ClippablePolygon;

/**
 * @brief G4Hype is a tube with hyperbolic profile; it describes an hyperbolic
 * volume with curved sides parallel to the Z axis. The solid has a specified
 * half-length along the Z axis, about which it is centered, and a given
 * minimum and maximum radii.
 */

class G4Hype : public G4VSolid
{
  public:

    /**
     * Constructs a hyperbolic tube, given its parameters.
     *  @param[in] pName The solid name.
     *  @param[in] newInnerRadius Inner radius.
     *  @param[in] newOuterRadius Outer radius.
     *  @param[in] newInnerStereo Inner stereo angle in radians.
     *  @param[in] newOuterStereo Outer stereo angle in radians.
     *  @param[in] newHalfLenZ Half length in Z.
     */
    G4Hype(const G4String& pName,
                 G4double newInnerRadius,
                 G4double newOuterRadius,
                 G4double newInnerStereo,
                 G4double newOuterStereo,
                 G4double newHalfLenZ);

    /**
     * Destructor.
     */
    ~G4Hype() override;

    /**
     * Accessors.
     */
    inline G4double GetInnerRadius () const;
    inline G4double GetOuterRadius () const;
    inline G4double GetZHalfLength () const;
    inline G4double GetInnerStereo () const;
    inline G4double GetOuterStereo () const;

    /**
     * Modifiers.
     */
    void SetInnerRadius (G4double newIRad);
    void SetOuterRadius (G4double newORad);
    void SetZHalfLength (G4double newHLZ);
    void SetInnerStereo (G4double newISte);
    void SetOuterStereo (G4double newOSte);

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
    EInside Inside(const G4ThreeVector& p) const override;
    G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;
    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const override;
    G4double DistanceToIn(const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;

    /**
     * Returns the type ID, "G4Hype" of the solid.
     */
    G4GeometryType  GetEntityType() const override;

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
    G4VisExtent GetExtent() const override;
    G4Polyhedron* CreatePolyhedron() const override;
    G4Polyhedron* GetPolyhedron() const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4Hype(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4Hype(const G4Hype& rhs);
    G4Hype& operator=(const G4Hype& rhs);

  private:

    /**
     * Tells whether we have an inner surface or not.
     */
    inline G4bool InnerSurfaceExists() const;

    /**
     * Returns the approximate isotropic distance to the hyperbolic surface.
     */
    G4double ApproxDistOutside( G4double pr, G4double pz,
                                G4double r0, G4double tanPhi ) const;
    G4double ApproxDistInside( G4double pr, G4double pz,
                               G4double r0, G4double tan2Phi ) const;

    /**
     * Returns the values of the hype radii at a given Z.
     */
    inline G4double HypeInnerRadius2(G4double zVal) const;
    inline G4double HypeOuterRadius2(G4double zVal) const;

    /**
     * Decides if and where a line intersects with a hyperbolic surface
     * (of infinite extent).
     *  @returns The number of intersections. If 0, the trajectory misses.
     */
    G4int IntersectHype( const G4ThreeVector& p, const G4ThreeVector& v,
                         G4double r2, G4double tan2Phi, G4double s[2] ) const;

  private:

    G4double innerRadius;
    G4double outerRadius;
    G4double halfLenZ;
    G4double innerStereo;
    G4double outerStereo;

    /** Precalculated parameters, squared quantities. */
    G4double tanInnerStereo;  // tan of Inner Stereo angle
    G4double tanOuterStereo;  // tan of Outer Stereo angle
    G4double tanInnerStereo2; // squared tan of Inner Stereo angle
    G4double tanOuterStereo2; // squared tan of Outer Stereo angle
    G4double innerRadius2;    // squared Inner Radius
    G4double outerRadius2;    // squared Outer Radius
    G4double endInnerRadius2; // squared endcap Inner Radius
    G4double endOuterRadius2; // squared endcap Outer Radius
    G4double endInnerRadius;  // endcap Inner Radius
    G4double endOuterRadius;  // endcap Outer Radius

    /** Used by DistanceToOut(). */
    enum ESide { outerFace, innerFace, leftCap, rightCap };

    G4double fCubicVolume = 0.0;
    G4double fSurfaceArea = 0.0;
    G4double fInnerSurfaceArea = 0.0;
    G4double fOuterSurfaceArea = 0.0;

    G4double fHalfTol;

    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;
};

#include "G4Hype.icc"

#endif  // defined(G4GEOM_USE_UHYPE) && defined(G4GEOM_USE_SYS_USOLIDS)

#endif // G4HYPE_HH
