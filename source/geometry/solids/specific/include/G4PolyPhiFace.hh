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
// G4PolyPhiFace
//
// Class description:
//
// Definition of a face that bounds a polycone or polyhedra when
// it has a phi opening:
//
//   G4PolyPhiFace( const G4ReduciblePolygon* rz,
//                        G4double phi,
//                        G4double deltaPhi,
//                        G4double phiOther )
//
// Specifically: a face that lies on a plane that passes through
// the z axis. It has boundaries that are straight lines of arbitrary
// length and direction, but with corners aways on the same side of
// the z axis.

// Author: David C. Williams (UCSC), 1998
// --------------------------------------------------------------------
#ifndef G4POLYPHIFACE_HH
#define G4POLYPHIFACE_HH 1

#include "G4VCSGface.hh"
#include "G4TwoVector.hh"

class G4ReduciblePolygon;

struct G4PolyPhiFaceVertex
{
  G4double x, y, r, z;   // position
  G4double rNorm, 
           zNorm;        // r/z normal
  G4ThreeVector norm3D;  // 3D normal

  // Needed for Triangulation Algorithm
  //
  G4bool ear;
  G4PolyPhiFaceVertex *next,*prev;
};

struct G4PolyPhiFaceEdge
{
  G4PolyPhiFaceEdge() = default;
  G4PolyPhiFaceVertex  *v0{nullptr}, *v1{nullptr};  // Corners
  G4double tr{.0}, tz{0.},        // Unit vector along edge
           length{0.};            // Length of edge
  G4ThreeVector norm3D;           // 3D edge normal vector
};

/**
 * @brief G4PolyPhiFace is a face that bounds a polycone or polyhedra when
 * it has a phi opening. Specifically, it is a face that lies on a plane that
 * passes through the Z axis, having boundaries that are straight lines of
 * arbitrary length and direction, but with corners aways on the same side of
 * the Z axis.
 */

class G4PolyPhiFace : public G4VCSGface
{

  public:

    /**
     * Constructor where points r,z should be supplied in clockwise order
     * in r,z.
     * For example:
     *                [1]---------[2]         ^ R
     *                 |           |          |
     *                 |           |          +--> z
     *                [0]---------[3]
     *  @param[in] rz Pointer to previous r,Z section.
     *  @param[in] phi Initial Phi starting angle.
     *  @param[in] deltaPhi Total Phi angle.
     *  @param[in] phiOther Phi angle of next section.
     */
    G4PolyPhiFace( const G4ReduciblePolygon* rz,
                         G4double phi, G4double deltaPhi, G4double phiOther );

    /**
     * Destructor. Removes edges and corners.
     */
    ~G4PolyPhiFace() override;

    /**
     * Copy constructor and assignment operator.
     */
    G4PolyPhiFace( const G4PolyPhiFace& source );
    G4PolyPhiFace& operator=( const G4PolyPhiFace& source );

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
                          const G4VoxelLimits &voxelLimit,
                          const G4AffineTransform& tranform,
                                G4SolidExtentList& extentList ) override;

    /**
     * Method invoked by the copy constructor or the assignment operator.
     * Its purpose is to return a pointer to a duplicate copy of the face.
     */
    inline G4VCSGface* Clone() override;

    /**
     * Returning an estimation of the face surface area, in internal units.
     */
    G4double SurfaceArea() override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4PolyPhiFace(__void__&);

    /**
     * Throws an exception if something is found inconsistent with the solid.
     * For debugging purposes only.
     */
    void Diagnose( G4VSolid* solid );

  private:

    /**
     * Calculates the surface area of a triangle. 
     * At the same time a random point in the triangle is given.
     */
    G4double SurfaceTriangle( const G4ThreeVector& p1, const G4ThreeVector& p2,
                              const G4ThreeVector& p3, G4ThreeVector* p4);

    /**
     * Auxiliary method for GetPointOnSurface().
     */
    G4ThreeVector GetPointOnFace() override;

    /**
     * Decides if the point in r,z is inside the edges of a face,
     * **but** do so consistently with other faces.
     */
    G4bool InsideEdgesExact( G4double r, G4double z, G4double normSign,
                             const G4ThreeVector& p, const G4ThreeVector& v );

    /**
     * Methods to decide if the point in r,z is inside the edges of a face.
     */
    G4bool InsideEdges( G4double r, G4double z );
    G4bool InsideEdges( G4double r, G4double z, G4double* distRZ2,
                        G4PolyPhiFaceVertex** base3Dnorm = nullptr,
                        G4ThreeVector** head3Dnorm = nullptr );

    /**
     * Decides precisely whether a trajectory passes to the left, right,
     * or exactly passes through the Z position of a vertex point in face.
     */
    inline G4double ExactZOrder( G4double z, 
                                 G4double qx, G4double qy, G4double qz, 
                           const G4ThreeVector& v, 
                                 G4double normSign,
                           const G4PolyPhiFaceVertex* vert ) const;

    /**
     * Copies parameters from other object; used in copy constructor and
     * assignment operator.
     */
    void CopyStuff( const G4PolyPhiFace& source );

    // Functions used for Triangulation in Case of generic Polygone.
    // The triangulation is used for GetPointOnFace()

    /**
     * Calculates of 2*Area of Triangle with Sign.
     */
    G4double Area2( const G4TwoVector& a, const G4TwoVector& b, const G4TwoVector& c);

    /**
     * Boolean functions for sign of Surface.
     */
    G4bool Left( const G4TwoVector& a, const G4TwoVector& b, const G4TwoVector& c );
    G4bool LeftOn( const G4TwoVector& a, const G4TwoVector& b, const G4TwoVector& c );
    G4bool Collinear( const G4TwoVector& a, const G4TwoVector& b, const G4TwoVector& c );

    /**
     * Boolean function for finding proper intersection of two
     * line segments (a,b) and (c,d).
     */
    G4bool IntersectProp( const G4TwoVector& a, const G4TwoVector& b,
                          const G4TwoVector& c, const G4TwoVector& d );

    /**
     * Boolean function for determining if point c is between a and b
     * where the three points (a,b,c) are on the same line.
     */
    G4bool Between( const G4TwoVector& a, const G4TwoVector& b, const G4TwoVector& c );

    /**
     * Boolean function for finding proper intersection or not
     * of two line segments (a,b) and (c,d).
     */
    G4bool Intersect( const G4TwoVector& a, const G4TwoVector& b,
                      const G4TwoVector& c, const G4TwoVector& d );

    /**
     * Boolean Diagonalie help to determine if diagonal s
     * of segment (a,b) is convex or reflex.
     */
    G4bool Diagonalie( G4PolyPhiFaceVertex* a, G4PolyPhiFaceVertex* b );

    /**
     * Boolean function for determining if b is inside the cone (a0,a,a1)
     * where a is the center of the cone.
     */
    G4bool InCone( G4PolyPhiFaceVertex *a, G4PolyPhiFaceVertex *b );

    /**
     * Boolean function for determining if Diagonal is possible
     * inside Polycone or PolyHedra.
     */
    G4bool Diagonal( G4PolyPhiFaceVertex* a, G4PolyPhiFaceVertex* b );

    /**
     * Initialisation for Triangulisation by ear tips.
     * For details see "Computational Geometry in C" by Joseph O'Rourke.
     */
    void EarInit();

    /**
     * Triangularisation by ear tips for Polycone or Polyhedra.
     * For details see "Computational Geometry in C" by Joseph O'Rourke.
     * NOTE: a copy of the shape is made and this copy is reordered in
     *       order to have a list of triangles. This list is used by the
     *       method GetPointOnFace().
     */
    void Triangulate();

  private:

    G4int numEdges = 0;  // Number of edges
    G4PolyPhiFaceEdge* edges = nullptr;       // The edges of the face
    G4PolyPhiFaceVertex* corners = nullptr;   // And the corners
    G4ThreeVector normal;        // Normal unit vector
    G4ThreeVector radial;        // Unit vector along radial direction
    G4ThreeVector surface;       // Point on surface
    G4ThreeVector surface_point; // Auxiliary point on surface used for
                                 // method GetPointOnFace() 
    G4double rMin, rMax, // Extent in r
             zMin, zMax; // Extent in z
    G4bool allBehind = false; // True if the polycone/polyhedra
                              // is behind the place of this face
    G4double kCarTolerance;      // Surface thickness
    G4double fSurfaceArea = 0.0; // Surface Area of PolyPhiFace 

     /** Auxiliary pointer to 'corners' used for triangulation.
         Copy structure, changing the structure of 'corners' (ear removal). */
    G4PolyPhiFaceVertex* triangles = nullptr;
};

#include "G4PolyPhiFace.icc"

#endif
