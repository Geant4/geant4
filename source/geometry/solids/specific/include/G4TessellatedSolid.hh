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
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * subject to DEFCON 705 IPR conditions.                            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4TessellatedSolid
//
// Class description:
//
// G4TessellatedSolid is a special Geant4 solid defined by a number of
// facets (UVFacet). It is important that the supplied facets shall form a
// fully enclose space which is the solid.
// At the moment only two types of facet can be used for the construction of
// a G4TessellatedSolid, i.e. the G4TriangularFacet and G4QuadrangularFacet.
//
// How to contruct a G4TessellatedSolid:
//
//    First declare a tessellated solid:
//
//      G4TessellatedSolid* solidTarget = new G4TessellatedSolid("Solid_name");
//
//    Define the facets which form the solid:
//
//      G4double targetSiz = 10*cm ;
//      G4TriangularFacet *facet1 = new
//      G4TriangularFacet (G4ThreeVector(-targetSize,-targetSize,        0.0),
//                         G4ThreeVector(+targetSize,-targetSize,        0.0),
//                         G4ThreeVector(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      G4TriangularFacet *facet2 = new
//      G4TriangularFacet (G4ThreeVector(+targetSize,-targetSize,        0.0),
//                         G4ThreeVector(+targetSize,+targetSize,        0.0),
//                         G4ThreeVector(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      G4TriangularFacet *facet3 = new
//      G4TriangularFacet (G4ThreeVector(+targetSize,+targetSize,        0.0),
//                         G4ThreeVector(-targetSize,+targetSize,        0.0),
//                         G4ThreeVector(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      G4TriangularFacet *facet4 = new
//      G4TriangularFacet (G4ThreeVector(-targetSize,+targetSize,        0.0),
//                         G4ThreeVector(-targetSize,-targetSize,        0.0),
//                         G4ThreeVector(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      G4QuadrangularFacet *facet5 = new
//      G4QuadrangularFacet (G4ThreeVector(-targetSize,-targetSize,      0.0),
//                           G4ThreeVector(-targetSize,+targetSize,      0.0),
//                           G4ThreeVector(+targetSize,+targetSize,      0.0),
//                           G4ThreeVector(+targetSize,-targetSize,      0.0),
//                           ABSOLUTE);
//
//    Then add the facets to the solid:
//
//      solidTarget->AddFacet((UVFacet*) facet1);
//      solidTarget->AddFacet((UVFacet*) facet2);
//      solidTarget->AddFacet((UVFacet*) facet3);
//      solidTarget->AddFacet((UVFacet*) facet4);
//      solidTarget->AddFacet((UVFacet*) facet5);
//
//    Finally declare the solid is complete:
//
//      solidTarget->SetSolidClosed(true);

// Author: P.R.Truscott (QinetiQ Ltd, UK), 31.10.2004 - Created.
//         M.Gayer (CERN), 12.10.2012 - New implementation with voxelization.
// --------------------------------------------------------------------
#ifndef G4TESSELLATEDSOLID_HH
#define G4TESSELLATEDSOLID_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UTESSELLATEDSOLID 1
#endif

#if defined(G4GEOM_USE_UTESSELLATEDSOLID)
  #define G4UTessellatedSolid G4TessellatedSolid
  #include "G4UTessellatedSolid.hh"
#else

#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "G4Types.hh"
#include "G4VSolid.hh"
#include "G4Voxelizer.hh"
#include "G4VFacet.hh"

struct G4VertexInfo
{
  G4int id;
  G4double mag2;
};

class G4VertexComparator
{
  public:

    G4bool operator() (const G4VertexInfo& l, const G4VertexInfo& r) const
    {
      return l.mag2 == r.mag2 ? l.id < r.id : l.mag2 < r.mag2;
    }
};

/**
 * @brief G4TessellatedSolid is a solid defined by a number of facets.
 * It is important that the supplied facets shall form a fully enclose space
 * which is the solid. The facets can be of two types, G4TriangularFacet and
 * G4QuadrangularFacet.
 */

class G4TessellatedSolid : public G4VSolid
{
  public:

    /**
     * Default Constructor.
     */
    G4TessellatedSolid ();

    /**
     * Constructor with solid's name.
     *  @param[in] name The name of the solid.
     */
    G4TessellatedSolid (const G4String& name);

    /**
     * Destructor. Clearing all allocated facets and data.
     */
    ~G4TessellatedSolid () override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4TessellatedSolid(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4TessellatedSolid (const G4TessellatedSolid& ts);
    G4TessellatedSolid& operator= (const G4TessellatedSolid& right);

    /**
     * Operator +=, allowing to add two tessellated solids together, so
     * that the solid on the left includes all of the facets in the solid
     * on the right. To note that copies of the facets are generated, rather
     * than using the original facet set of the solid on the right.
     */
    G4TessellatedSolid& operator+= (const G4TessellatedSolid& right);

    /**
     * Methods for adding or retrieving a facet given an index.
     */
    G4bool AddFacet (G4VFacet* aFacet);
    inline G4VFacet* GetFacet (G4int i) const;

    /**
     * Accessors.
     */
    G4int GetNumberOfFacets () const;
    G4int GetFacetIndex (const G4ThreeVector& p) const;
    G4double GetMinXExtent () const;
    G4double GetMaxXExtent () const;
    G4double GetMinYExtent () const;
    G4double GetMaxYExtent () const;
    G4double GetMinZExtent () const;
    G4double GetMaxZExtent () const;


    /**
     * Concrete implementations of the expected query interfaces for
     * solids, as defined in the base class G4VSolid.
     */
    EInside Inside (const G4ThreeVector& p) const override;
    G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;
    G4double DistanceToIn(const G4ThreeVector& p,
                                  const G4ThreeVector& v)const override;
    G4double DistanceToIn(const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p,
                                   const G4ThreeVector& v,
                                   const G4bool calcNorm,
                                         G4bool* validNorm,
                                         G4ThreeVector* norm) const override;

    /**
     * Returns the outwards pointing unit normal of the shape for the
     * surface closest to the point at offset 'p'.
     *  @param[in] p The point coordinates.
     *  @param[out] n The returned normal vector.
     *  @returns false if not a valid normal.
     */
    virtual G4bool Normal (const G4ThreeVector& p, G4ThreeVector& n) const;

    /**
     * Returns the the safety distance from outside the solid at a point 'p'.
     *  @param[in] p The point coordinates.
     *  @param[in] aAccurate Accuracy flag, if false quickly computes and
     *              returns the distance to the voxels bounding-box.
     *  @returns The safety distance.
     */
    virtual G4double SafetyFromOutside(const G4ThreeVector& p,
                                             G4bool aAccurate = false) const;

    /**
     * Returns the the safety distance from inside the solid at a point 'p'.
     *  @param[in] p The point coordinates.
     *  @param[in] aAccurate Not used.
     *  @returns The safety distance.
     */
    virtual G4double SafetyFromInside (const G4ThreeVector& p,
                                             G4bool aAccurate = false) const;

    /**
     * Returns the type ID, "G4TessellatedSolid" of the solid.
     */
    G4GeometryType GetEntityType () const override;

    /**
     * Returns true as the solid has only planar faces.
     */
    G4bool IsFaceted () const override;

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo(std::ostream& os) const override;

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
     * Returning an estimation of the solid volume (capacity) and
     * surface area, in internal units.
     */
    G4double GetSurfaceArea() override;
    G4double GetCubicVolume() override;

    /**
     * Modifier and accessor to close/finalise the solid.
     */
    void SetSolidClosed (const G4bool t);
    G4bool GetSolidClosed () const;

    /**
     * Checks the structure of the solid.
     *  @returns A value, sum of the following defect indicators, if any
     *           (0 means no defects):
     *           1 - cubic volume is negative, wrong orientation of facets;
     *           2 - some facets have wrong orientation;
     *           4 - holes in the surface.
     */
    G4int CheckStructure() const;

    /**
     * Allowing to tune the maximum number of voxels to use for optimisation.
     */
    inline void SetMaxVoxels(G4int max);

    /**
     * Returns the voxels structure.
     */
    inline G4Voxelizer& GetVoxels();

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
     * Computes the bounding limits of the solid.
     *  @param[out] pMin The minimum bounding limit point.
     *  @param[out] pMax The maximum bounding limit point.
     */
    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    G4Polyhedron* CreatePolyhedron() const override;
    G4Polyhedron* GetPolyhedron() const override;
    void DescribeYourselfTo(G4VGraphicsScene& scene) const override;
    G4VisExtent GetExtent() const override;

    /**
     * Loggers reporting the total allocated memory.
     */
    G4int AllocatedMemoryWithoutVoxels();
    G4int AllocatedMemory();
    void DisplayAllocatedMemory();

  private:

    /**
     * Initialisation/reset of data, used in constructors and operators.
     */
    void Initialize();

    /**
     * Resetting/copying data, used in constructors and operators.
     */
    void DeleteObjects ();
    void CopyObjects (const G4TessellatedSolid& s);

    /**
     * Internal methods used for computing distances with or without voxels.
     */
    G4double DistanceToOutNoVoxels(const G4ThreeVector& p,
                                   const G4ThreeVector& v,
                                         G4ThreeVector& aNormalVector,
                                         G4bool&        aConvex,
                                         G4double aPstep = kInfinity) const;
    G4double DistanceToInCandidates(const std::vector<G4int>& candidates,
                                    const G4ThreeVector& aPoint,
                                    const G4ThreeVector& aDirection) const;
    void DistanceToOutCandidates(const std::vector<G4int>& candidates,
                                 const G4ThreeVector& aPoint,
                                 const G4ThreeVector& direction,
                                       G4double& minDist,
                                       G4ThreeVector& minNormal,
                                       G4int& minCandidate) const;
    G4double DistanceToInNoVoxels(const G4ThreeVector& p,
                                  const G4ThreeVector& v,
                                        G4double aPstep = kInfinity) const;
    G4double DistanceToInCore(const G4ThreeVector &p, const G4ThreeVector& v,
                                    G4double aPstep = kInfinity) const;
    G4double DistanceToOutCore(const G4ThreeVector& p, const G4ThreeVector& v,
                                     G4ThreeVector& aNormalVector,
                                     G4bool& aConvex,
                                     G4double aPstep = kInfinity) const;

    /**
     * Finds those facets that have surface planes that bound the volume.
     * To note that this is going to reject concave surfaces as being extreme.
     */
    void SetExtremeFacets();

    /**
     * Internal methods used for checking if a point 'p' is inside the solid
     * in presence or not of voxels.
     */
    EInside InsideNoVoxels (const G4ThreeVector& p) const;
    EInside InsideVoxels(const G4ThreeVector& p) const;

    /**
     * Performs the voxelisation of the shape, building the optimisation
     * structure, according to the specified parameters.
     */
    void Voxelize();

    /**
     * Creates a list of vertices with an additional sorted list, where all
     * the items are sorted by magnitude of vertices vector.
     */
    void CreateVertexList();

    /**
     * Utilities for preparation of voxels indeces. Used in Voxelize() function.
     */
    void PrecalculateInsides();
    G4int SetAllUsingStack(const std::vector<G4int>& voxel,
                           const std::vector<G4int>& max,
                                 G4bool status, G4SurfBits& checked);

    /**
     * Utility to compare sorted voxels.
     */
    static G4bool CompareSortedVoxel(const std::pair<G4int, G4double>& l,
                                     const std::pair<G4int, G4double>& r);


    /**
      * Prepares a set of predefined random vectors, used to generate rays
      * from a user-defined point. Used in Inside() function to determine
      * whether the point is inside or outside of the tessellated solid.
      * All vectors should be unit vectors.
      */
    void SetRandomVectors();

    /**
     * Computes the minimum distance of a point 'p' from a 'facet'.
     */
    G4double MinDistanceFacet(const G4ThreeVector& p, G4bool simple,
                                    G4VFacet* &facet) const;

    /**
     * Computes if a point 'p' is outside or not of the computed extent,
     * given a 'tolerance'. Used internally in Inside() functions.
     *  @returns true if the point is within the extent.
     */
    inline G4bool OutsideOfExtent(const G4ThreeVector& p,
                                        G4double tolerance = 0.0) const;

  protected:

    G4double kCarToleranceHalf;

  private:

    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;

    std::vector<G4VFacet*> fFacets;
    std::set<G4VFacet*> fExtremeFacets; // Does all other facets lie on
                                        // or behind this surface?

    G4GeometryType fGeometryType;
    G4double       fCubicVolume = 0.0;
    G4double       fSurfaceArea = 0.0;

    std::vector<G4ThreeVector> fVertexList;

    std::set<G4VertexInfo,G4VertexComparator> fFacetList;

    G4ThreeVector fMinExtent, fMaxExtent;

    G4bool fSolidClosed = false;

    std::vector<G4ThreeVector> fRandir;

    G4int fMaxTries;

    G4Voxelizer fVoxels;  // Pointer to the voxelized solid

    G4SurfBits fInsides;
};

///////////////////////////////////////////////////////////////////////////////
// Inline Methods
///////////////////////////////////////////////////////////////////////////////

inline G4VFacet* G4TessellatedSolid::GetFacet (G4int i) const
{
  return fFacets[i];
}

inline void G4TessellatedSolid::SetMaxVoxels(G4int max)
{
  fVoxels.SetMaxVoxels(max);
}

inline G4Voxelizer& G4TessellatedSolid::GetVoxels()
{
  return fVoxels;
}

inline G4bool G4TessellatedSolid::OutsideOfExtent(const G4ThreeVector& p,
                                                  G4double tolerance) const
{
  return ( p.x() < fMinExtent.x() - tolerance
        || p.x() > fMaxExtent.x() + tolerance
        || p.y() < fMinExtent.y() - tolerance
        || p.y() > fMaxExtent.y() + tolerance
        || p.z() < fMinExtent.z() - tolerance
        || p.z() > fMaxExtent.z() + tolerance);
}

#endif

#endif
