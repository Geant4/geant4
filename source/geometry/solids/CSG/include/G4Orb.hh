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
// G4Orb
//
// Class description:
//
// A G4Orb represents a full sphere; it is a simple case of G4Sphere.

// Author: Vladimir Grichine (CERN), 20.08.2003 - Created
//         Evgueni Tcherniaev (CERN), 08.08.2017 - Revised
// --------------------------------------------------------------------
#ifndef G4ORB_HH
#define G4ORB_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UORB 1
#endif

#if defined(G4GEOM_USE_UORB)
  #define G4UOrb G4Orb
  #include "G4UOrb.hh"
#else

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4CSGSolid.hh"
#include "G4Polyhedron.hh"

/**
 * @brief G4Orb represents a full sphere.
 */

class G4Orb : public G4CSGSolid
{
  public:

    /**
     * Constructs a full sphere, given a name and its radius.
     *  @param[in] pName The name of the solid.
     *  @param[in] pRmax Outer radius.
     */
    G4Orb(const G4String& pName, G4double pRmax);

    /**
     * Default destructor.
     */
    ~G4Orb() override = default;

    /**
     * Accessors and modifiers.
     */
    inline G4double GetRadius() const;
    inline G4double GetRadialTolerance() const;
    inline void SetRadius(G4double newRmax);

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
     * Returns the type ID, "G4Orb" of the solid.
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
    void DescribeYourselfTo (G4VGraphicsScene& scene) const override;
    G4VisExtent GetExtent() const override;
    G4Polyhedron* CreatePolyhedron() const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4Orb(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4Orb(const G4Orb& rhs) = default;
    G4Orb& operator=(const G4Orb& rhs);

  protected:

    /**
     * Checks radius and initialises data members. Used in constructor.
     */
    void Initialize();

  private:

    G4double fRmax = 0.0;
    G4double halfRmaxTol = 0.0;
    G4double sqrRmaxPlusTol = 0.0, sqrRmaxMinusTol = 0.0;
};

#include "G4Orb.icc"

#endif

#endif
