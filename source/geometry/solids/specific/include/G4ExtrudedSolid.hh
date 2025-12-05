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
// G4ExtrudedSolid
//
// Class description:
//
// G4ExtrudedSolid is a solid which represents the extrusion of an arbitrary
// polygon with fixed outline in the defined Z sections.
// The z-sides of the solid are the scaled versions of the same polygon.
// The solid is implemented as a specification of G4TessellatedSolid.
//
// Parameters in the constructor:
// const G4String& pName             - solid name
// std::vector<G4TwoVector> polygon  - the vertices of the outlined polygon
//                                     defined in clockwise or anti-clockwise order
// std::vector<ZSection>             - the z-sections defined by
//                                     z position, offset and scale
//                                     in increasing z-position order
//
// Parameters in the special constructor (for solid with 2 z-sections:
// G4double halfZ                    - the solid half length in Z
// G4TwoVector off1                  - offset of the side in -halfZ
// G4double scale1                   - scale of the side in -halfZ
// G4TwoVector off2                  - offset of the side in +halfZ
// G4double scale2                   - scale of the side in +halfZ

// Author: Ivana Hrivnacova (IPN, Orsay), 09.02.2007 - First implementation
// --------------------------------------------------------------------
#ifndef G4EXTRUDEDSOLID_HH
#define G4EXTRUDEDSOLID_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UEXTRUDEDSOLID 1
#endif

#if defined(G4GEOM_USE_UEXTRUDEDSOLID)
  #define G4UExtrudedSolid G4ExtrudedSolid
  #include "G4UExtrudedSolid.hh"
#else

#include <vector>

#include "G4TwoVector.hh"
#include "G4TessellatedSolid.hh"

/**
 * @brief G4ExtrudedSolid is a is a solid which represents the extrusion
 * of an arbitrary polygon with fixed outline in the defined Z sections.
 * The z-sides of the solid are the scaled versions of the same polygon.
 * The solid is implemented as a specification of a G4TessellatedSolid.
 */

class G4ExtrudedSolid : public G4TessellatedSolid
{
  public:

    /**
     * Structure defining a Z section composing the solid.
     */
    struct ZSection
    {
      ZSection() : fZ(0.), fOffset(0.,0.), fScale(1.) {}
      ZSection(G4double z, const G4TwoVector& offset, G4double scale)
        : fZ(z), fOffset(offset), fScale(scale) {}

      G4double    fZ;
      G4TwoVector fOffset;
      G4double    fScale;
    };

    /**
     * General constructor for an extruded polygon, through contour and polyline.
     *  @param[in] pName The solid name.
     *  @param[in] polygon The 2D polygonal contour, i.e. the vertices of the
     *             outlined polygon defined in clock-wise order.
     *  @param[in] zsections The 3D polyline with scale factors, i.e. the
     *             Z-sections defined by Z position in increasing order.
     */
    G4ExtrudedSolid( const G4String&                 pName,
                     const std::vector<G4TwoVector>& polygon,
                     const std::vector<ZSection>&    zsections);

    /**
     * Special constructor for an extruded polygon with 2 Z-sections.
     *  @param[in] pName The solid name.
     *  @param[in] polygon The 2D polygonal contour, i.e. the vertices of the
     *             outlined polygon defined in clock-wise order.
     *  @param[in] halfZ Half length in Z, i.e. the distance from the origin
     *             to the sections.
     *  @param[in] off1 (X, Y) position of the first polygon in -halfZ.
     *  @param[in] scale1 Scale factor at -halfZ.
     *  @param[in] off2 (X, Y) position of the second polygon in +halfZ.
     *  @param[in] scale2 Scale factor at +halfZ.
     */
    G4ExtrudedSolid( const G4String&                 pName,
                     const std::vector<G4TwoVector>& polygon,
                           G4double                  halfZ,
                     const G4TwoVector& off1 = G4TwoVector(0.,0.),
                           G4double scale1 = 1.,
                     const G4TwoVector& off2 = G4TwoVector(0.,0.),
                           G4double scale2 = 1. );

    /**
     * Default Destructor.
     */
    ~G4ExtrudedSolid() override = default;

    /**
     * Accessors.
     */
    inline G4int       GetNofVertices() const;
    inline G4TwoVector GetVertex(G4int index) const;
    inline std::vector<G4TwoVector> GetPolygon() const;
    inline G4int       GetNofZSections() const;
    inline ZSection    GetZSection(G4int index) const;
    inline std::vector<ZSection> GetZSections() const;

    /**
     * Concrete implementations of the expected query interfaces for
     * solids, as defined in the base class G4VSolid.
     */
    EInside  Inside(const G4ThreeVector& p) const override;
    G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;
    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const override;
    G4double DistanceToIn(const G4ThreeVector& p ) const override;
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
                                 G4double& pMin, G4double& pMax) const override;

    /**
     * Returns the type ID, "G4ExtrudedSolid" of the solid.
     */
    G4GeometryType GetEntityType () const override;

    /**
     * Returns true as the solid has only planar faces.
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
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4ExtrudedSolid(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4ExtrudedSolid(const G4ExtrudedSolid& rhs) = default;
    G4ExtrudedSolid& operator=(const G4ExtrudedSolid& rhs);

  private:

    /**
     * Algorithm for SurfaceNormal() following the original
     * specification for points not on the surface.
     */
    G4ThreeVector ApproxSurfaceNormal(const G4ThreeVector& p) const;

    /**
     * Computes parameters for point projections p(z)
     * to the polygon scale & offset.
     */
    void ComputeProjectionParameters();

    /**
     * Computes the lateral planes: a*x + b*y + c*z + d = 0.
     */
    void ComputeLateralPlanes();

    /**
     * Returns if point 'p' is within the polygon.
     */
    inline G4bool PointInPolygon(const G4ThreeVector& p) const;

    /**
     * Returns the square distance of point 'p' from the polygon.
     */
    inline G4double DistanceToPolygonSqr(const G4ThreeVector& p) const;

    /**
     * Returns the vertex coordinates, given the indeces for the
     * polygons and Z sections.
     *  @param[in] iz Index for the Z section.
     *  @param[in] ind Index for the polygon.
     *  @returns The shifted and scaled coordinates of the vertex.
     */
    G4ThreeVector GetVertex(G4int iz, G4int ind) const;

    /**
     * Returns the projected point of 'p' in the polygon scale.
     */
    G4TwoVector ProjectPoint(const G4ThreeVector& point) const;

    /**
     * Returns true if 'p' is on the line through 'l1', 'l2'.
     */
    G4bool IsSameLine(const G4TwoVector& p,
                      const G4TwoVector& l1,
                      const G4TwoVector& l2) const;
    /**
     * Returns true if 'p' is on the line through 'l1', 'l2'
     * and lies between 'l1' and 'l2'.
     */
    G4bool IsSameLineSegment(const G4TwoVector& p,
                             const G4TwoVector& l1,
                             const G4TwoVector& l2) const;
    /**
     * Returns true if 'p1' and 'p2' are on the same side of the line
     * through 'l1', 'l2'.
     */
    G4bool IsSameSide(const G4TwoVector& p1,
                      const G4TwoVector& p2,
                      const G4TwoVector& l1,
                      const G4TwoVector& l2) const;
    /**
     * Returns true if 'p' is inside of triangle abc or on its edges.
     */
    G4bool IsPointInside(const G4TwoVector& a,
                         const G4TwoVector& b,
                         const G4TwoVector& c,
                         const G4TwoVector& p) const;
    /**
     * Returns the angle of the vertex in 'p0'.
     */
    G4double GetAngle(const G4TwoVector& p0,
                      const G4TwoVector& pa,
                      const G4TwoVector& pb) const;

    /**
     * Returns a pointer to a triangular facet from the polygon points
     * given by indices forming the down side ( the normal goes in -z).
     */
    G4VFacet* MakeDownFacet(G4int ind1, G4int ind2, G4int ind3) const;

    /**
     * Returns a pointer to a triangular facet from the polygon points
     * given by indices forming the upper side ( z>0 ).
     */
    G4VFacet* MakeUpFacet(G4int ind1, G4int ind2, G4int ind3) const;

    /**
     * Decomposes polygonal sides in triangular facets.
     *  @returns false if failing to define a facet.
     */
    G4bool AddGeneralPolygonFacets();

    /**
     * Generates the tessellated structure of the solid creating the
     * triangular or quadrangular facets from the vertices.
     *  @returns false if failing to define a facet.
     */
    G4bool MakeFacets();

  private:

    std::size_t    fNv;
    std::size_t    fNz;
    std::vector<G4TwoVector> fPolygon;
    std::vector<ZSection>    fZSections;
    std::vector< std::vector<G4int> > fTriangles;
    G4bool         fIsConvex = false;
    G4GeometryType fGeometryType;

    G4int fSolidType = 0;
    struct plane { G4double a,b,c,d; }; // a*x + b*y + c*z + d = 0
    std::vector<plane> fPlanes;
    struct line { G4double k,m; };      // x = k*y + m;
    std::vector<line> fLines;
    std::vector<G4double> fLengths;     // edge lengths

    std::vector<G4double>      fKScales;
    std::vector<G4double>      fScale0s;
    std::vector<G4TwoVector>   fKOffsets;
    std::vector<G4TwoVector>   fOffset0s;
};

#include "G4ExtrudedSolid.icc"

#endif

#endif
