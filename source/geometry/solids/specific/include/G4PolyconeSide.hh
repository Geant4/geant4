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
// G4PolyconeSide
//
// Class description:
//
// Class implementing a face that represents one conical side of a polycone:
//
//   G4PolyconeSide( const G4PolyconeSideRZ* prevRZ,
//                   const G4PolyconeSideRZ* tail,
//                   const G4PolyconeSideRZ* head,
//                   const G4PolyconeSideRZ* nextRZ,
//                         G4double phiStart, G4double deltaPhi, 
//                         G4bool phiIsOpen, G4bool isAllBehind=false )
//
// Values for r1,z1 and r2,z2 should be specified in clockwise order in (r,z).

// Author: David C. Williams (UCSC), 1998 - Created
// --------------------------------------------------------------------
#ifndef G4POLYCONESIDE_HH
#define G4POLYCONESIDE_HH

#include "G4VCSGface.hh"

class G4IntersectingCone;

struct G4PolyconeSideRZ
{
  G4double r, z;  // start of vector
};

// ----------------------------------------------------------------------------
// MT-specific utility code 

#include "G4GeomSplitter.hh"

// The class G4PlSideData is introduced to encapsulate the
// fields of the class G4PolyconeSide that may not be read-only.
//
class G4PlSideData
{
  public:

    void initialize()
    {
      fPhix = 0.; fPhiy = 0.; fPhiz = 0.; fPhik = 0.;
    }

    G4double fPhix=0., fPhiy=0., fPhiz=0., fPhik=0.;   // Cached values for phi
};

// The type G4PlSideManager is introduced to 
// encapsulate the methods used by both the master thread and 
// worker threads to allocate memory space for the fields encapsulated
// by the class G4PlSideData.
//
using G4PlSideManager = G4GeomSplitter<G4PlSideData>;

//
// ----------------------------------------------------------------------------

/**
 * @brief G4PolyconeSide is a utility class implementing a face that
 * represents one conical side of a polycone.
 */

class G4PolyconeSide : public G4VCSGface
{
  public:

    /**
     * Constructor for the conical side of a polycone.
     *  @param[in] prevRZ Pointer to previous r,Z section.
     *  @param[in] tail Pointer to r,Z tail of section.
     *  @param[in] head Pointer to r,Z head of section.
     *  @param[in] nextRZ Pointer to next r,Z section.
     *  @param[in] phiStart Initial Phi starting angle.
     *  @param[in] deltaPhi Total Phi angle.
     *  @param[in] phiIsOpen Flag indicating if it is a Phi section.
     *  @param[in] isAllBehind Indicating if entire surface is behind normal.
     */
    G4PolyconeSide( const G4PolyconeSideRZ* prevRZ,
                    const G4PolyconeSideRZ* tail,
                    const G4PolyconeSideRZ* head,
                    const G4PolyconeSideRZ* nextRZ,
                          G4double phiStart, G4double deltaPhi, 
                          G4bool phiIsOpen, G4bool isAllBehind = false );

    /**
     * Destructor.
     */
    ~G4PolyconeSide() override;
  
    /**
     * Copy constructor and assignment operator.
     */
    G4PolyconeSide( const G4PolyconeSide& source );
    G4PolyconeSide& operator=( const G4PolyconeSide& source );
  
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
    G4bool Intersect(const G4ThreeVector& p, const G4ThreeVector& v,  
                           G4bool outgoing, G4double surfTolerance,
                           G4double& distance, G4double &distFromSurface,
                           G4ThreeVector& normal, G4bool& isAllBehind) override;

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
    EInside Inside( const G4ThreeVector& p, G4double tolerance, 
                          G4double* bestDistance ) override;
  
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
                          const G4VoxelLimits& voxelLimit,
                          const G4AffineTransform& tranform,
                                G4SolidExtentList& extentList ) override;

    /**
     * Method invoked by the copy constructor or the assignment operator.
     * Its purpose is to return a pointer to a duplicate copy of the face.
     */
    inline G4VCSGface* Clone() override { return new G4PolyconeSide( *this ); }

    /**
     * Returning an estimation of the face surface area, in internal units.
     */
    G4double SurfaceArea() override;

    /**
     * Returns a random point located and uniformly distributed on the face.
     */
    G4ThreeVector GetPointOnFace() override;
  
    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4PolyconeSide(__void__&);

    /**
     * Returns the instance ID.
     */
    inline G4int GetInstanceID() const  { return instanceID; }

    /**
     * Returns the private data instance manager.
     */
    static const G4PlSideManager& GetSubInstanceManager();

  private:

    /**
     * Calculates the distance of a point from the conical surface, including
     * the effect of any phi segmentation.
     *  @param[in] p The point to check.
     *  @param[in] opposite If true, check the opposite hemisphere (see below).
     *  @param[out] distOutside Additional distance outside the edges of surface.
     *  @param[out] rzNorm If negative, the point is inside.
     *  @returns The distance from the conical plane, if extrapolated beyond
     *           edges, signed by whether the point is in inside or outside
     *           the shape.
     */
    G4double DistanceAway( const G4ThreeVector& p, G4bool opposite,
                                 G4double& distOutside,
                                 G4double* rzNorm = nullptr );
      
    /**
     * Special version of DistanceAway() for Inside. Opposite parameter is not
     * used, instead use sign of rx for choosing the side.
     *  @param[in] p The point to check.
     *  @param[out] distOutside Additional distance outside the edges of surface.
     *  @param[out] edgeRZnorm If negative, the point is inside.
     *  @returns The distance from the conical plane.
     */
    G4double DistanceAway( const G4ThreeVector& p, G4double& distOutside,
                                 G4double* edgeRZnorm );

    /**
     * Decides if a point is on a cone and returns the 'normal' if it is.
     *  @returns true if the point is on the cone.
     */
    G4bool PointOnCone( const G4ThreeVector& hit, G4double normSign,
                        const G4ThreeVector& p,
                        const G4ThreeVector& v, G4ThreeVector& normal );

    /**
     * Copies parameters from other object; used in copy constructor and
     * assignment operator.
     */
    void CopyStuff( const G4PolyconeSide& source );
  
    /**
     * Decides the point at which two 2-dimensional lines intersect.
     * It is assumed that the lines are *not* parallel.
     */
    static void FindLineIntersect( G4double x1, G4double y1,
                                   G4double tx1, G4double ty1,
                                   G4double x2, G4double y2,
                                   G4double tx2, G4double ty2,
                                   G4double& x, G4double& y );

    /**
     * Calculates Phi for a given 3-vector (point 'p'), if not already cached
     * for the same point, in the attempt to avoid consecutive computation of
     * the same quantity.
     */
    G4double GetPhi( const G4ThreeVector& p );

  private:

    G4double r[2], z[2]; // r, z parameters, in specified order
    G4double startPhi,   // Start phi (0 to 2pi), if phiIsOpen
             deltaPhi;   // Delta phi (0 to 2pi), if phiIsOpen
    G4bool phiIsOpen = false; // True if there is a phi slice
    G4bool allBehind = false; // True if the entire solid is "behind" this face
  
    G4IntersectingCone* cone = nullptr;  // Our intersecting utility class
  
    G4double rNorm, zNorm;  // Normal to surface in r,z space
    G4double rS, zS;        // Unit vector along surface in r,z space
    G4double length;        // Length of face in r,z space
    G4double prevRS,
             prevZS;        // Unit vector along previous polyconeSide
    G4double nextRS,
             nextZS;        // Unit vector along next polyconeSide
  
    G4double rNormEdge[2],
             zNormEdge[2];  // Normal to edges

    G4int ncorners = 0;
    G4ThreeVector* corners = nullptr; // The coordinates of the corners
                                      // (if phiIsOpen)

    G4double kCarTolerance;       // Geometrical surface thickness
    G4double fSurfaceArea = 0.0;  // Used for surface calculation 

    G4int instanceID;
      // This field is used as instance ID.
    G4GEOM_DLL static G4PlSideManager subInstanceManager;
      // This field helps to use the class G4PlSideManager introduced above.
};

#endif
