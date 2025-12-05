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
// G4Torus
//
// Class description:
//
// A torus or torus segment with curved sides parallel to the z-axis.
// The torus has a specified swept radius about which it is centered,
// and a given minimum and maximum radius. A minimum radius of 0
// signifies a filled torus.
// The torus segment is specified by starting and delta angles for phi,
// with 0 being the +x axis, PI/2 the +y axis. A delta angle of 2PI
// signifies a complete, unsegmented torus/cylinder.
//
// Member functions:
//
//   As inherited from G4CSGSolid+
//
//     G4Torus(const G4String      &pName
//             G4double      pRmin
//             G4double      pRmax
//             G4double      pRtor
//             G4double      pSPhi
//             G4double      pDPhi )
//
//     - Construct a torus with the given name and dimensions.
//       The angles are provided is radians. pRtor >= pRmax
//
// Member Data:
//
//   fRmin  Inside radius
//   fRmax  Outside radius
//   fRtor  swept radius of torus
//
//   fSPhi  The starting phi angle in radians,
//          adjusted such that fSPhi+fDPhi<=2PI, fSPhi>-2PI
//
//   fDPhi  Delta angle of the segment in radians
//
// You could find very often in G4Torus functions values like 'pt' or
// 'it'. These are the distances from p or i G4ThreeVector points in the
// plane (Z axis points p or i) to fRtor point in XY plane. This value is
// similar to rho for G4Tubs and is used for definiton of the point
// relative to fRmin and fRmax, i.e. for solution of inside/outside
// problems

// Author: V.Grichine (CERN), 30.10.1996 - First version
//         E.Medernach (CERN), 31.08.2000 - Migrated to numeric solutions
// --------------------------------------------------------------------
#ifndef G4TORUS_HH
#define G4TORUS_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UTORUS 1
#endif

#if (defined(G4GEOM_USE_UTORUS) && defined(G4GEOM_USE_SYS_USOLIDS))
  #define G4UTorus G4Torus
  #include "G4UTorus.hh"
#else

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4CSGSolid.hh"

/**
 * @brief G4Torus represents a torus or torus segment with curved sides
 * parallel to the z-axis. The torus has a specified swept radius about which
 * it is centered, and a given minimum and maximum radius. A minimum radius
 * of 0 signifies a filled torus.
 * The torus segment is specified by starting and delta angles for phi,
 * with 0 being the +x axis, PI/2 the +y axis. A delta angle of 2PI
 * signifies a complete, unsegmented torus/cylinder.
 */

class G4Torus : public G4CSGSolid
{

  public:

    /**
     * Constructs a torus or torus segment with the given name and dimensions.
     *  @param[in] pName The name of the solid.
     *  @param[in] pRmin Inner radius.
     *  @param[in] pRmax Outer radius.
     *  @param[in] pRtor Swept radius of torus.
     *  @param[in] pSPhi Starting Phi angle in radians
     *             adjusted such that fSPhi+fDPhi<=2PI, fSPhi>-2PI.
     *  @param[in] pDPhi Delta angle of the segment in radians.
     */
    G4Torus(const G4String& pName,
                  G4double pRmin,
                  G4double pRmax,
                  G4double pRtor,
                  G4double pSPhi,
                  G4double pDPhi);

    /**
     * Default destructor.
     */
    ~G4Torus() override = default;

    /**
     * Accessors.
     */
    inline G4double GetRmin() const;
    inline G4double GetRmax() const;
    inline G4double GetRtor() const;
    inline G4double GetSPhi() const;
    inline G4double GetDPhi() const;
    inline G4double GetSinStartPhi () const;
    inline G4double GetCosStartPhi () const;
    inline G4double GetSinEndPhi   () const;
    inline G4double GetCosEndPhi   () const;

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
    void ComputeDimensions(      G4VPVParameterisation* p,
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
    G4double DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;

    /**
     * Returns the type ID, "G4Torus" of the solid.
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
    G4Polyhedron* CreatePolyhedron () const override;

    /**
     * Checks and sets all the parameters given in input. Used in constructor.
     */
    void SetAllParameters(G4double pRmin, G4double pRmax, G4double pRtor,
                          G4double pSPhi, G4double pDPhi);

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4Torus(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4Torus(const G4Torus& rhs) = default;
    G4Torus& operator=(const G4Torus& rhs);

  private:

    /**
     * Calculates the real roots to the torus surface, using the
     * G4JTPolynomialSolver class. Returns negative solutions as well.
     */
    void TorusRootsJT(const G4ThreeVector& p,
                      const G4ThreeVector& v,
                            G4double r,
                            std::vector<G4double>& roots) const ;

    /**
     * Interface method for DistanceToIn() and DistanceToOut().
     * Calls TorusRootsJT() using the Jenkins-Traub algorithm for real
     * polynomial root finding.
     *  @returns The smalles possible distance to the surface.
     */
    G4double SolveNumericJT(const G4ThreeVector& p,
                            const G4ThreeVector& v,
                                  G4double r,
                                  G4bool IsDistanceToIn) const;

    /**
     * Algorithm for SurfaceNormal() following the original specification
     * for points not on the surface.
     */
    G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p) const;

  private:

    /** The cached parameters, ensured within range. */
    G4double fRmin, fRmax, fRtor, fSPhi, fDPhi;

    /** Radial and angular tolerances. */
    G4double fRminTolerance, fRmaxTolerance, kRadTolerance, kAngTolerance;

    /** Cached half tolerance values. */
    G4double halfCarTolerance, halfAngTolerance;
};

#include "G4Torus.icc"

#endif  // defined(G4GEOM_USE_UTORUS) && defined(G4GEOM_USE_SYS_USOLIDS)


#endif // G4TORUS_HH
