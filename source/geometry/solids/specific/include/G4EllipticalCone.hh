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
// G4EllipticalCone
//
// Class description:
//
// G4EllipticalCone is a full cone with elliptical base which can be cut in Z.
//
// Member Data:
//
//      xSemiAxis       semi-axis, x, without dimentions
//      ySemiAxis       semi-axis, y, without dimentions
//      zheight         height, z
//      zTopCut         upper cut plane level, z
//
// The height in Z corresponds to where the elliptical cone hits the
// Z-axis if it had no Z cut. Also the cone is centered at zero having a
// base at zTopCut and another at -zTopCut. The semi-major axes at the Z=0
// plane are given by xSemiAxis*zheight and ySemiAxis*zheight so that the
// curved surface of our cone satisfies the equation:
//
// ***************************************************************************
// *                                                                         *
// *           (x/xSemiAxis)^2 + (y/ySemiAxis)^2 = (zheight - z)^2           *
// *                                                                         *
// ***************************************************************************
//
// In case you want to construct G4EllipticalCone from:
//   1. halflength in Z = zTopCut
//   2. Dx and Dy =  halflength of ellipse axis  at  z = -zTopCut
//   3. dx and dy =  halflength of ellipse axis  at  z =  zTopCut
//      ! Attention :  dx/dy=Dx/Dy
//
// You need to find xSemiAxis,ySemiAxis and zheight:
//
//  xSemiAxis = (Dx-dx)/(2*zTopCut)
//  ySemiAxis = (Dy-dy)/(2*zTopCut)
//    zheight = (Dx+dx)/(2*xSemiAxis)

// Author: Dionysios Anninos (CERN), 08.09.2005 - Created
//         Lukas Lindroos, Tatiana Nikitina (CERN), 20.08.2007 - Revised
//         Evgueni Tcherniaev (CERN), 20.07.2017 - New revision
// --------------------------------------------------------------------
#ifndef G4ELLIPTICALCONE_HH
#define G4ELLIPTICALCONE_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UELLIPTICALCONE 1
#endif

#if (defined(G4GEOM_USE_UELLIPTICALCONE) && defined(G4GEOM_USE_SYS_USOLIDS))
  #define G4UEllipticalCone G4EllipticalCone
  #include "G4UEllipticalCone.hh"
#else

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4VSolid.hh"
#include "G4Polyhedron.hh"

/**
 * @brief G4EllipticalCone is a full cone with elliptical base which
 * can be cut in Z. The height in Z corresponds to where the elliptical
 * cone hits the Z-axis if it had no Z cut.
 */

class G4EllipticalCone : public G4VSolid
{
  public:

    /**
     * Constructs an elliptical cone, with cut in Z.
     *  @param[in] name The solid name.
     *  @param[in] pxSemiAxis Scalar value, defining the scaling along X-axis.
     *  @param[in] pySemiAxis Scalar value, defining the scaling along Y-axis.
     *  @param[in] zMax The Z-coordinate at the apex.
     *  @param[in] pzTopCut Upper cut plane level.
     */
    G4EllipticalCone(const G4String& pName,
                           G4double  pxSemiAxis,
                           G4double  pySemiAxis,
                           G4double  zMax,
                           G4double  pzTopCut);

    /**
     * Destructor.
     */
    ~G4EllipticalCone() override;

    /**
     * Accessors.
     */
    inline G4double GetSemiAxisMin () const;
    inline G4double GetSemiAxisMax () const;
    inline G4double GetSemiAxisX () const;
    inline G4double GetSemiAxisY () const;
    inline G4double GetZMax() const;
    inline G4double GetZTopCut() const;

    /**
     * Modifiers.
     */
    void SetSemiAxis (G4double x, G4double y, G4double z);
    void SetZCut (G4double newzTopCut);

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
    G4double DistanceToOut(const G4ThreeVector& p,
                           const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;

    /**
     * Returns the type ID, "G4EllipticalCone" of the solid.
     */
    G4GeometryType GetEntityType() const override;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Returns a random point located and uniformly distributed on the
     * surface of the solid.
     */
    G4ThreeVector GetPointOnSurface() const override;

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo(std::ostream& os) const override;

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    G4Polyhedron* GetPolyhedron() const override;
    void DescribeYourselfTo(G4VGraphicsScene& scene) const override;
    G4VisExtent GetExtent() const override;
    G4Polyhedron* CreatePolyhedron() const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4EllipticalCone(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4EllipticalCone(const G4EllipticalCone& rhs);
    G4EllipticalCone& operator=(const G4EllipticalCone& rhs);

  protected:

    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;

  private:

    /**
     * Algorithm for SurfaceNormal() following the original
     * specification for points not on the surface.
     */
    G4ThreeVector ApproxSurfaceNormal(const G4ThreeVector& p) const;

  private:

    G4double halfCarTol;
    G4double fCubicVolume = 0.0;
    G4double fSurfaceArea = 0.0;
    G4double fMinZBaseArea = 0.0;
    G4double fMaxZBaseArea = 0.0;
    G4double xSemiAxis, ySemiAxis, zheight, zTopCut;
    G4double cosAxisMin, invXX, invYY;
};

#include "G4EllipticalCone.icc"

#endif  // defined(G4GEOM_USE_UELLIPTICALCONE) && defined(G4GEOM_USE_SYS_USOLIDS)

#endif // G4ELLIPTICALCONE_HH
