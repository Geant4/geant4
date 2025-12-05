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
// G4EllipticalTube
//
// Class description:
//
// A tube with elliptical cross section:
//
//   G4EllipticalTube( const G4String& name,
//                           G4double  Dx,
//                           G4double  Dy,
//                           G4double  Dz )
//
// The equation of the lateral surface is: (x/dx)^2 + (y/dy)^2 = 1

// Author: David C. Williams (UCSC), 29.03.2000 - First implementation
//         Evgueni Tcherniaev (CERN), 23.12.2019 - Revised
// --------------------------------------------------------------------
#ifndef G4ELLIPTICALTUBE_HH
#define G4ELLIPTICALTUBE_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UELLIPTICALTUBE 1
#endif

#if (defined(G4GEOM_USE_UELLIPTICALTUBE) && defined(G4GEOM_USE_SYS_USOLIDS))
  #define G4UEllipticalTube G4EllipticalTube
  #include "G4UEllipticalTube.hh"
#else

#include "G4VSolid.hh"
#include "G4Polyhedron.hh"

/**
 * @brief G4EllipticalTube is a tube with elliptical cross section.
 * The equation of the lateral surface is: (x/dx)^2 + (y/dy)^2 = 1.
 */

class G4EllipticalTube : public G4VSolid
{
  public:

    /**
     * Constructs an elliptical tube, given its parameters.
     *  @param[in] name The solid name.
     *  @param[in] Dx Half length of axis along X.
     *  @param[in] Dy Half length of axis along Y.
     *  @param[in] Dz Half length in Z.
     */
    G4EllipticalTube( const G4String& name,
                            G4double Dx,
                            G4double Dy,
                            G4double Dz );

    /**
     * Destructor.
     */
    ~G4EllipticalTube() override;

    /**
     * Accessors.
     */
    inline G4double GetDx() const;
    inline G4double GetDy() const;
    inline G4double GetDz() const;

    /**
     * Modifiers.
     */
    inline void SetDx( G4double Dx );
    inline void SetDy( G4double Dy );
    inline void SetDz( G4double Dz );

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
     * Returns the type ID, "G4EllipticalTube" of the solid.
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
    G4Polyhedron* CreatePolyhedron() const override;
    G4Polyhedron* GetPolyhedron() const override;
    void DescribeYourselfTo( G4VGraphicsScene& scene ) const override;
    G4VisExtent GetExtent() const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4EllipticalTube(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4EllipticalTube(const G4EllipticalTube& rhs);
    G4EllipticalTube& operator=(const G4EllipticalTube& rhs);

  private:

    /**
     * Checks parameters and sets pre-calculated values.
     */
    void CheckParameters();

    /**
     * Algorithm for SurfaceNormal() following the original
     * specification for points not on the surface.
     */
    G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p ) const;

    /**
     * Calculates the surface area and caches it.
     */
    G4double GetCachedSurfaceArea() const;

  private:

    G4double halfTolerance;

    G4double fDx; // semi-axis in X
    G4double fDy; // semi-axis in Y
    G4double fDz; // half length in Z

    G4double fCubicVolume = 0.0; // volume
    G4double fSurfaceArea = 0.0; // surface area

    /** Cached pre-calculated values. */
    G4double fRsph;    // R of bounding sphere
    G4double fDDx;     // Dx squared
    G4double fDDy;     // Dy squared
    G4double fSx;      // X scale factor
    G4double fSy;      // Y scale factor
    G4double fR;       // resulting Radius, after scaling elipse to circle
    G4double fQ1;      // distance approximation : dist = Q1*(x^2 + y^2) - Q2
    G4double fQ2;      // distance approximation : dist = Q1*(x^2 + y^2) - Q2
    G4double fScratch; // half length of scratching segment squared

    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;
};

#include "G4EllipticalTube.icc"

#endif  // defined(G4GEOM_USE_UELLIPTICALTUBE) && defined(G4GEOM_USE_SYS_USOLIDS)

#endif // G4ELLIPTICALTUBE_HH
