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
// G4GenericTrap
//
// Class description:
//
// G4GenericTrap is a solid which represents an arbitrary trapezoid with
// up to 8 vertices standing on two parallel planes perpendicular to Z axis.
//
// Parameters in the constructor:
// - name               - solid name
// - halfZ              - the solid half length in Z
// - vertices           - the (x,y) coordinates of vertices:
//                        o first four points: vertices[i], i<4
//                          are the vertices sitting on the -halfZ plane;
//                        o last four points: vertices[i], i>=4
//                          are the vertices sitting on the +halfZ plane.
//
// The order of defining the vertices of the solid is the following:
//      - point 0 is connected with points 1,3,4
//      - point 1 is connected with points 0,2,5
//      - point 2 is connected with points 1,3,6
//      - point 3 is connected with points 0,2,7
//      - point 4 is connected with points 0,5,7
//      - point 5 is connected with points 1,4,6
//      - point 6 is connected with points 2,5,7
//      - point 7 is connected with points 3,4,6
// Points can be identical in order to create shapes with less than 8 vertices.
// Adapted from Arb8 implementation in Root/TGeo.

// Authors: T.Nikitina (CERN) & I.Hrivnacova (IPN, Orsay), 27.05.2010 - Created
//          Evgueni Tcherniaev (CERN), 27.05.2024 - Complete revision, speed up
// -------------------------------------------------------------------
#ifndef G4GENERICTRAP_HH
#define G4GENERICTRAP_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UGENERICTRAP 1
#endif

#if defined(G4GEOM_USE_UGENERICTRAP)
  #define G4UGenericTrap G4GenericTrap
  #include "G4UGenericTrap.hh"
#else

#include <vector>

#include "globals.hh"
#include "G4TwoVector.hh"
#include "G4VSolid.hh"

/**
 * @brief G4GenericTrap is a solid which represents an arbitrary trapezoid with
 * up to 8 vertices standing on two parallel planes perpendicular to the Z axis.
 * Points can be identical in order to create shapes with less than 8 vertices.
 */

class G4GenericTrap : public G4VSolid
{
  public:

    /**
     * Constructs an generic trapezoid, given its vertices.
     *  @param[in] name The solid name.
     *  @param[in] halfZ Half length in Z.
     *  @param[in] vertices The (x,y) coordinates of the vertices.
     */
    G4GenericTrap(const G4String& name, G4double halfZ,
                  const std::vector<G4TwoVector>& vertices);

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4GenericTrap(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4GenericTrap(const G4GenericTrap& rhs);
    G4GenericTrap& operator=(const G4GenericTrap& rhs);

    /**
     * Default Destructor.
     */
    ~G4GenericTrap() override = default;

    /**
     * Accessors and modifiers.
     */
    inline G4double    GetZHalfLength() const;
    inline G4int       GetNofVertices() const;
    inline G4TwoVector GetVertex(G4int index) const;
    inline const std::vector<G4TwoVector>& GetVertices() const;
    inline G4double    GetTwistAngle(G4int index) const;
    inline G4bool      IsTwisted() const;
    inline G4int       GetVisSubdivisions() const;
    inline void        SetVisSubdivisions(G4int subdiv);

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
     * Returns the type ID, "G4GenericTrap" of the solid.
     */
    G4GeometryType GetEntityType() const override;

    /**
     * Returns true if the solid has only planar faces; false if twisted.
     */
    G4bool IsFaceted () const override;

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
     * Returns a random point located and uniformly distributed on the
     * surface of the solid.
     */
    G4ThreeVector GetPointOnSurface() const override ;

    /**
     * Returning an estimation of the solid volume (capacity) and
     * surface area, in internal units.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo(G4VGraphicsScene& scene) const override;
    G4VisExtent GetExtent() const override;
    G4Polyhedron* CreatePolyhedron() const override;
    G4Polyhedron* GetPolyhedron() const override;

  private:

    /**
     * Algorithm for SurfaceNormal() following the original
     * specification for points not on the surface.
     */
    G4ThreeVector ApproxSurfaceNormal(const G4ThreeVector& p) const;

    /**
     * Checks the parameters of the solid and issues exception if leading
     * to an invalid construct.
     */
    void CheckParameters(G4double halfZ, const std::vector<G4TwoVector>& vertices);

    /**
     * Computes surface equations and twist angles of lateral faces.
     */
    void ComputeLateralSurfaces();

    /**
     * Sets the bounding box.
     */
    void ComputeBoundingBox();

    /**
     * Sets the max length of a scratch.
     */
    void ComputeScratchLength();

    /**
     * Computes the lateral face area, given the face index.
     * Used for random sampling of points on surface.
     */
    G4double GetLateralFaceArea(G4int iface) const;

    /**
     * Logger methods for issuing warnings.
     */
    void WarningSignA(const G4String& method, const G4String& icase, G4double A,
                      const G4ThreeVector& p, const G4ThreeVector& v) const;
    void WarningSignB(const G4String& method, const G4String& icase, G4double f, G4double B,
                      const G4ThreeVector& p, const G4ThreeVector& v) const;
    void WarningDistanceToIn(G4int k, const G4ThreeVector& p, const G4ThreeVector& v,
                             G4double tmin, G4double tmax,
                             const G4double ttin[2], const G4double ttout[2]) const;
    void WarningDistanceToOut(const G4ThreeVector& p,
                              const G4ThreeVector& v,
                              G4double tout) const;

  private:

    struct G4GenericTrapPlane // Ax + By + Cz + D = 0
    {
      G4double A = 0.;
      G4double B = 0.;
      G4double C = 0.;
      G4double D = 0.;
    };
    struct G4GenericTrapSurface // Axz + Byz + Czz + Dx + Ey + Fz + G = 0
    {
      G4double A = 0.;
      G4double B = 0.;
      G4double C = 0.;
      G4double D = 0.;
      G4double E = 0.;
      G4double F = 0.;
      G4double G = 0.;
    };

    // Data members
    G4double                 halfTolerance = 0.;
    G4double                 fScratch = 0.;
    G4double                 fDz = 0.;
    std::vector<G4TwoVector> fVertices = {0.,0.,0.,0.,0.,0.,0.,0.};
    G4TwoVector              fDelta[4];
    G4bool                   fIsTwisted = false;
    G4double                 fTwist[5] = {0.};
    G4ThreeVector            fMinBBox{0.};
    G4ThreeVector            fMaxBBox{0.};
    G4int                    fVisSubdivisions = 0;
    G4GenericTrapPlane       fPlane[8];
    G4GenericTrapSurface     fSurf[4];
    G4double                 f4k[4] = {0.}; // Lipschitz constants * 4
    mutable G4double         fArea[4] = {0.};
    mutable G4bool           fRebuildPolyhedron = false;
    mutable G4Polyhedron*    fpPolyhedron = nullptr;

    // Surface and Volume
    G4double                 fSurfaceArea = 0.;
    G4double                 fCubicVolume = 0.;
};

#include "G4GenericTrap.icc"

#endif // defined(G4GEOM_USE_UGENERICTRAP)

#endif // G4GENERICTRAP_HH
