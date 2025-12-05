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
// G4PolyhedraSide
//
// Class description:
//
// Class implementing a face that represents one segmented side of a polyhedra:
//
//   G4PolyhedraSide( const G4PolyhedraSideRZ* prevRZ,
//                    const G4PolyhedraSideRZ* tail,
//                    const G4PolyhedraSideRZ* head,
//                    const G4PolyhedraSideRZ* nextRZ,
//                          G4int numSide,
//                          G4double phiStart, G4double phiTotal, 
//                          G4bool phiIsOpen, G4bool isAllBehind=false )
//
// Values for r1,z1 and r2,z2 should be specified in clockwise order in (r,z).

// Author: David C. Williams (UCSC), 1998 - Created
// --------------------------------------------------------------------
#ifndef G4POLYHEDRASIDE_HH
#define G4POLYHEDRASIDE_HH

#include "G4VCSGface.hh"

class G4IntersectingCone;

struct G4PolyhedraSideRZ
{
  G4double r, z;  // start of vector
};

// ----------------------------------------------------------------------------
// MT-specific utility code 

#include "G4GeomSplitter.hh"

// The class G4PhSideData is introduced to encapsulate the
// fields of the class G4PolyhedraSide that may not be read-only.
//
class G4PhSideData
{
  public:

    void initialize()
    {
      fPhix = 0.; fPhiy = 0.; fPhiz = 0.; fPhik = 0.;
    }

    G4double fPhix=0., fPhiy=0., fPhiz=0., fPhik=0.;   // Cached values for phi
};

// The type G4PhSideManager is introduced to encapsulate the methods used
// by both the master thread and worker threads to allocate memory space
// for the fields encapsulated by the class G4PhSideData.
//
using G4PhSideManager = G4GeomSplitter<G4PhSideData>;

//
// ----------------------------------------------------------------------------

/**
 * @brief G4PolyhedraSide is a utility class implementing a face that
 * represents one segmented side of a polyhedra.
 */

class G4PolyhedraSide : public G4VCSGface
{

  public:

    /**
     * Constructor for the segmented side of a polyhedra.
     *  @param[in] prevRZ Pointer to previous r,Z section.
     *  @param[in] tail Pointer to r,Z tail of section.
     *  @param[in] head Pointer to r,Z head of section.
     *  @param[in] nextRZ Pointer to next r,Z section.
     *  @param[in] numSide The number od sides.
     *  @param[in] phiStart Initial Phi starting angle.
     *  @param[in] phiTotal Total Phi angle.
     *  @param[in] phiIsOpen Flag indicating if it is a Phi section.
     *  @param[in] isAllBehind Indicating if entire surface is behind normal.
     */
    G4PolyhedraSide( const G4PolyhedraSideRZ* prevRZ,
                     const G4PolyhedraSideRZ* tail,
                     const G4PolyhedraSideRZ* head,
                     const G4PolyhedraSideRZ* nextRZ,
                           G4int numSide,
                           G4double phiStart, G4double phiTotal, 
                           G4bool phiIsOpen, G4bool isAllBehind = false );

    /**
     * Destructor.
     */
    ~G4PolyhedraSide() override;
  
    /**
     * Copy constructor and assignment operator.
     */
    G4PolyhedraSide( const G4PolyhedraSide& source );
    G4PolyhedraSide& operator=( const G4PolyhedraSide& source );
  
    /**
     * Determines the distance along a line to the face.
     *  @param[in] p Position.
     *  @param[in] v Direction (assumed to be a unit vector).
     *  @param[in] outgoing Flag true, to consider only inside surfaces;
     *             false, to consider only outside surfaces.
     *  @param[in] surfTolerance Minimum distance from the surface.
     *  @param[out] distance Distance to intersection.
     *  @param[out] distFromSurface Distance from surface (along surface normal),
     *              < 0 if the point is in front of the surface.
     *  @param[out] normal Normal of surface at intersection point.
     *  @param[out] allBehind Flag, true, if entire surface is behind normal.
     *  @returns true if there is an intersection, false otherwise.
     */
    G4bool Intersect( const G4ThreeVector& p, const G4ThreeVector& v,  
                            G4bool outgoing, G4double surfTolerance,
                            G4double& distance, G4double& distFromSurface,
                            G4ThreeVector& normal, G4bool& allBehind ) override;

    /**
     * Determines the distance of a point from either the inside or outside
     * surfaces of the face.
     *  @param[in] p Position.
     *  @param[in] outgoing Flag, true, to consider only inside surfaces
     *             or false, to consider only outside surfaces.
     *  @returns The distance to the closest surface satisfying requirements
     *           or kInfinity if no such surface exists.
     */
    G4double Distance( const G4ThreeVector& p, G4bool outgoing ) override;
  
    /**
     * Determines whether a point is inside, outside, or on the surface of
     * the face.
     *  @param[in] p Position.
     *  @param[in] tolerance Tolerance defining the bounds of the "kSurface",
     *             nominally equal to kCarTolerance/2.
     *  @param[out] bestDistance Distance to the closest surface (in or out).
     *  @returns kInside if the point is closest to the inside surface;
     *           kOutside if the point is closest to the outside surface;
     *           kSurface if the point is withing tolerance of the surface.
     */
    EInside Inside( const G4ThreeVector &p, G4double tolerance, 
                          G4double *bestDistance ) override;
  
    /**
     * Returns the normal of surface closest to the point.
     *  @param[in] p Position.
     *  @param[out] bestDistance Distance to the closest surface (in or out).
     *  @returns The normal of the surface nearest the point.
     */
    G4ThreeVector Normal( const G4ThreeVector& p,
                                G4double* bestDistance ) override;

    /**
     * Returns the face extent along the axis.
     *  @param[in] axis Unit vector defining the direction.
     *  @returns The largest point along the given axis of the face's extent.
     */
    G4double Extent( const G4ThreeVector axis ) override;
  
    /**
     * Calculates the extent of the face for the voxel navigator.
     *  @param[in] axis The axis in which to check the shapes 3D extent against.
     *  @param[in] voxelLimit Limits along x, y, and/or z axes.
     *  @param[in] tranform A coordinate transformation on which to apply to
     *             the shape before testing.
     *  @param[out] extentList The list of (voxel) extents along the axis.
     */
    void CalculateExtent( const EAxis axis, 
                          const G4VoxelLimits &voxelLimit,
                          const G4AffineTransform& tranform,
                                G4SolidExtentList& extentList ) override;

    /**
     * Method invoked by the copy constructor or the assignment operator.
     * Its purpose is to return a pointer to a duplicate copy of the face.
     */
    inline G4VCSGface* Clone() override { return new G4PolyhedraSide( *this ); }

    /**
     * Returning an estimation of the face surface area, in internal units.
     */
    G4double SurfaceArea() override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4PolyhedraSide(__void__&);

    /**
     * Returns the instance ID.
     */
    inline G4int GetInstanceID() const  { return instanceID; }

    /**
     * Returns the private data instance manager.
     */
    static const G4PhSideManager& GetSubInstanceManager();

  private:

    /**
     * Internal data structures.
     */
    struct sG4PolyhedraSideVec;         // Secret recipe for allowing
    friend struct sG4PolyhedraSideVec;  // protected nested structures

    using G4PolyhedraSideEdge = struct sG4PolyhedraSideEdge
    {
      G4ThreeVector  normal;       // Unit normal to this edge
      G4ThreeVector  corner[2];    // The two corners of this phi edge
      G4ThreeVector  cornNorm[2];  // The normals of these corners
    };
  
    using G4PolyhedraSideVec = struct sG4PolyhedraSideVec
    {
      G4ThreeVector  normal,   // Normal (point out of the shape)
                     center,   // Point in center of side
                     surfPhi,  // Unit vector on surface pointing along phi
                     surfRZ;   // Unit vector on surface pointing along R/Z
      G4PolyhedraSideEdge* edges[2];  // The phi boundary edges to this side 
                                      //     [0]=low phi [1]=high phi
      G4ThreeVector edgeNorm[2];      // RZ edge normals [i] at {r[i],z[i]}
    };

    /**
     * Calculates the surface area of a triangle. 
     * At the same time a random point in the triangle is given.
     */
    G4double SurfaceTriangle( const G4ThreeVector& p1,
                              const G4ThreeVector& p2,
                              const G4ThreeVector& p3,
                                    G4ThreeVector* p4 );

    /**
     * Returns a random point located and uniformly distributed on the face.
     */
    G4ThreeVector GetPointOnFace() override;  

    /**
     * Auxiliary method for GetPointOnSurface().
     */
    G4ThreeVector GetPointOnPlane( const G4ThreeVector& p0, const G4ThreeVector& p1, 
                                   const G4ThreeVector& p2, const G4ThreeVector& p3,
                                   G4double* Area );

    /**
     * Decides if a line correctly intersects one side plane of a segment.
     * It is assumed that the correct side has been chosen, and thus only 
     * the Z bounds (of the entire segment) are checked.
     *  @param[in] p The point to check.
     *  @param[in] v The direction.
     *  @param[in] vec Description record of the side plane.
     *  @param[in] normSign Sign (+/- 1) to apply to normal.
     *  @param[in] surfTolerance Surface tolerance (generally > 0, see below).
     *  @param[out] distance Distance along v to intersection.
     *  @param[out] distFromSurface Distance from surface normal.
     *  @returns true is the line is intersecting, false otherwise.
     */
    G4bool IntersectSidePlane( const G4ThreeVector& p, const G4ThreeVector& v,
                               const G4PolyhedraSideVec& vec,
                                     G4double normSign, 
                                     G4double surfTolerance,
                                     G4double& distance,
                                     G4double& distFromSurface );

    /**
     * Calculates which phi segments a line intersects in three dimensions.
     * No check is made as to whether the intersections are within the Z bounds
     * of the segment.
     */
    G4int LineHitsSegments( const G4ThreeVector& p,
                            const G4ThreeVector& v,
                                  G4int* i1, G4int* i2 );

    /**
     * Decides which phi segment an angle belongs to, counting from zero.
     * A returned value of -1 indicates that the phi value is outside the
     * shape (only possible if phiTotal < 360 degrees).
     */
    G4int PhiSegment( G4double phi );

    /**
     * Decides which phi segment is closest in phi to the point.
     * The result is the same as PhiSegment() if there is no phi opening.
     */
    G4int ClosestPhiSegment( G4double phi );
  
    /**
     * Calculates Phi for a given 3-vector (point 'p'), if not already cached
     * for the same point, in the attempt to avoid consecutive computation of
     * the same quantity.
     */
    G4double GetPhi( const G4ThreeVector& p );

    /**
     * Computes the total distance from the side.
     *  @param[in] p The point to check.
     *  @param[in] vec The vector set of this side.
     *  @param[out] normDist The returned distance normal to the side or edge,
     *              as appropriate, signed.
     *  @returns The total distance from the side.
     */
    G4double DistanceToOneSide( const G4ThreeVector& p,
                                const G4PolyhedraSideVec& vec,
                                      G4double* normDist );

    /**
     * Calculates the distance of a point from the segmented surface, including
     * the effect of any phi segmentation.
     * Adds distance from side edges, if necessary, to the total distance,
     * and updates 'normDist' appropriate depending on edge normals.
     *  @param[in] p The point to check.
     *  @param[in] opposite If true, check the opposite hemisphere (see below).
     *  @param[out] normDist The returned normal distance.
     *  @returns The distance from the segmented plane.
     */
    G4double DistanceAway( const G4ThreeVector& p,
                           const G4PolyhedraSideVec& vec,
                                 G4double* normDist );
             
    /**
     * Copies parameters from other object; used in copy constructor and
     * assignment operator.
     */
    void CopyStuff( const G4PolyhedraSide& source );

  private:

    G4int numSide = 0;    // Number sides
    G4double r[2], z[2];  // r, z parameters, in specified order
    G4double startPhi,    // Start phi (0 to 2pi), if phiIsOpen
             deltaPhi,    // Delta phi (0 to 2pi), if phiIsOpen
             endPhi;      // End phi (>startPhi), if phiIsOpen
    G4bool phiIsOpen = false; // True if there is a phi slice
    G4bool allBehind = false; // True if the entire solid is "behind" this face
  
    G4IntersectingCone* cone = nullptr; // Our intersecting cone
  
    G4PolyhedraSideVec* vecs = nullptr; // Vector set for each facet of our face
    G4PolyhedraSideEdge* edges = nullptr; // The edges belong to vecs
    G4double    lenRZ,      // RZ length of each side
                lenPhi[2];  // Phi dimensions of each side
    G4double    edgeNorm;   // Normal in RZ/Phi space to each side

    G4double kCarTolerance;      // Geometrical surface thickness
    G4double fSurfaceArea = 0.0; // Surface Area 

    G4int instanceID;
      // This field is used as instance ID.
    G4GEOM_DLL static G4PhSideManager subInstanceManager;
      // This field helps to use the class G4PhSideManager introduced above.
};

#endif
