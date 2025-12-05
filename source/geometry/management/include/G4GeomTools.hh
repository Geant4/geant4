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
// G4GeomTools
//
// Class description:
//
// A collection of utilities which can be helpful for a wide range
// of geometry-related tasks

// Author: Evgueni Tcherniaev (CERN), 10.10.2016
// --------------------------------------------------------------------
#ifndef G4GEOMTOOLS_HH
#define G4GEOMTOOLS_HH

#include <vector>
#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"

using G4TwoVectorList = std::vector<G4TwoVector>;
using G4ThreeVectorList = std::vector<G4ThreeVector>;

/**
 * @brief G4GeomTools is a collecting utilities which can be helpful
 * for a wide range of geometry-related tasks.
 */

class G4GeomTools
{
  public:

    // ==================================================================
    //   2D Utilities
    // ------------------------------------------------------------------

    /**
     * Functions to calculate the area of 2D triangle, returned value is
     * positive if the vertices of the triangle are given in anticlockwise
     * order, otherwise it is negative.
     */
    static G4double TriangleArea(G4double Ax, G4double Ay,
                                 G4double Bx, G4double By,
                                 G4double Cx, G4double Cy);
    static G4double TriangleArea(const G4TwoVector& A,
                                 const G4TwoVector& B,
                                 const G4TwoVector& C);

    /**
     * Calculates the area of a 2D quadrilateral, returned value is positive if
     * the vertices of the quadrilateral are given in anticlockwise order,
     * otherwise it is negative.
     */
    static G4double QuadArea(const G4TwoVector& A,
                             const G4TwoVector& B,
                             const G4TwoVector& C,
                             const G4TwoVector& D);

    /**
     * Calculates the area of a 2D polygon, returned value is positive if
     * the vertices of the polygon are defined in anticlockwise order,
     * otherwise it is negative.
     */
    static G4double PolygonArea(const G4TwoVectorList& polygon);

    /**
     * Decides if a point (Px,Py) is inside the triangle (Ax,Ay)(Bx,By)(Cx,Cy).
     */
    static G4bool PointInTriangle(G4double Px, G4double Py,
                                  G4double Ax, G4double Ay,
                                  G4double Bx, G4double By,
                                  G4double Cx, G4double Cy);

    /**
     * Decides if a point P is inside the triangle ABC.
     */
    static G4bool PointInTriangle(const G4TwoVector& P,
                                  const G4TwoVector& A,
                                  const G4TwoVector& B,
                                  const G4TwoVector& C);

    /**
     * Decides if a point P is inside the 'Polygon'.
     */
    static G4bool PointInPolygon(const G4TwoVector& P,
                                 const G4TwoVectorList& Polygon);

    /**
     * Decides if a 2D 'polygon' is convex, i.e. if all internal angles are
     * less than pi.
     */
    static G4bool IsConvex(const G4TwoVectorList& polygon);

    /**
     * Simple implementation of "ear clipping" algorithm for triangulation
     * of a simple contour/polygon, it places results in a std::vector as
     * triplets of vertices. If triangulation is successful the function
     * returns true, otherwise false.
     */
    static G4bool TriangulatePolygon(const G4TwoVectorList& polygon,
                                           std::vector<G4int>& result);

    /**
     * Same using the function above and returning as 'result' a list
     * of triangles.
     */
    static G4bool TriangulatePolygon(const G4TwoVectorList& polygon,
                                           G4TwoVectorList& result);

    /**
     * Removes collinear and coincident points from a 2D 'polygon'.
     * Indices of removed points are available in 'iout'.
     * Allows to specify a 'tolerance'.
     */
    static void RemoveRedundantVertices(G4TwoVectorList& polygon,
                                        std::vector<G4int>& iout,
                                        G4double tolerance = 0.0);

    /**
     * Calculates the bounding rectangle of a disk sector. It returns false
     * if the input parameters do not meet the following criteria:
     *   rmin   >= 0
     *   rmax   >  rmin + kCarTolerance
     *   delPhi >  0 + kCarTolerance.
     */
    static G4bool DiskExtent(G4double rmin, G4double rmax,
                             G4double startPhi, G4double delPhi,
                             G4TwoVector& pmin, G4TwoVector& pmax);

    /**
     * Calculates the bounding rectangle of a disk sector.
     * Faster version without check of parameters.
     */
    static void DiskExtent(G4double rmin, G4double rmax,
                           G4double sinPhiStart, G4double cosPhiStart,
                           G4double sinPhiEnd, G4double cosPhiEnd,
                           G4TwoVector& pmin, G4TwoVector& pmax);

    /**
     * Computes the circumference (perimeter) of an ellipse.
     */
    static G4double EllipsePerimeter(G4double a,
                                     G4double b);

    /**
     * Computes the lateral surface area of an elliptic cone.
     */
    static G4double EllipticConeLateralArea(G4double a,
                                            G4double b,
                                            G4double h);

    // ==================================================================
    //   3D Utilities
    // ------------------------------------------------------------------

    /**
     * Finds the normal to the plane of a 3D triangle ABC;
     * the length of the normal is equal to the area of the triangle.
     */
    static G4ThreeVector TriangleAreaNormal(const G4ThreeVector& A,
                                            const G4ThreeVector& B,
                                            const G4ThreeVector& C);

    /**
     * Finds the normal to the plane of a 3D quadrilateral ABCD;
     * the length of the normal is equal to the area of the quadrilateral.
     */
    static G4ThreeVector QuadAreaNormal(const G4ThreeVector& A,
                                        const G4ThreeVector& B,
                                        const G4ThreeVector& C,
                                        const G4ThreeVector& D);

    /**
     * Finds the normal to the plane of a 3D polygon; the length of the
     * normal is equal to the area of the polygon.
     */
    static G4ThreeVector PolygonAreaNormal(const G4ThreeVectorList& polygon);

    /**
     * Calculates the distance between a point 'P' and line segment AB in 3D.
     */
    static G4double DistancePointSegment(const G4ThreeVector& P,
                                         const G4ThreeVector& A,
                                         const G4ThreeVector& B);

    /**
     * Finds a point on a 3D line segment AB closest to point 'P'.
     */
    static G4ThreeVector ClosestPointOnSegment(const G4ThreeVector& P,
                                               const G4ThreeVector& A,
                                               const G4ThreeVector& B);

    /**
     * Finds a point on a 3D triangle ABC closest to point 'P'.
     */
    static G4ThreeVector ClosestPointOnTriangle(const G4ThreeVector& P,
                                                const G4ThreeVector& A,
                                                const G4ThreeVector& B,
                                                const G4ThreeVector& C);

    /**
     * Calculates the bounding box of a spherical sector,
     *  @returns false if input parameters do not meet the following criteria:
     *   rmin       >= 0
     *   rmax       >  rmin + kCarTolerance
     *   startTheta >= 0 && <= pi;
     *   delTheta   >  0 + kCarTolerance
     *   delPhi     >  0 + kCarTolerance.
     */
    static G4bool SphereExtent(G4double rmin, G4double rmax,
                               G4double startTheta, G4double delTheta,
                               G4double startPhi, G4double delPhi,
                               G4ThreeVector& pmin, G4ThreeVector& pmax);

    /**
     * Calculates the hyperbolic surface stereo. Stereo is a half angle at the
     * intersection point of the two lines in the tangent plane cross-section.
     */
    static G4double HypeStereo(G4double r0, // radius at z = 0
                               G4double r,  // radius at z = h
                               G4double h);

    /**
     * Finds the XY-coordinates of the corners of a generic trap that bounds
     * specified twisted tube.
     *  @param[in] twistAng The twist angle.
     *  @param[in] endInnerRad The inner radius at z = halfZ.
     *  @param[in] endOuterRad The outer radius at z = halfZ.
     *  @param[in] dPhi Delta phi.
     *  @param[in] vertices The corners of the generic trap.
     */
    static void TwistedTubeBoundingTrap(G4double twistAng,
                                        G4double endInnerRad,
                                        G4double endOuterRad,
                                        G4double dPhi,
                                        G4TwoVectorList& vertices);

    /**
     * Calculate surface area of the hyperboloid between 'zmin' and 'zmax'.
     *  @param[in] dphi Delta phi.
     *  @param[in] r0 The radius at z = 0.
     *  @param[in] tanstereo The tangent of the stereo angle.
     *  @param[in] zmin Minimum Z.
     *  @param[in] zmax Maximum Z.
     */
    static G4double HyperboloidSurfaceArea(G4double dphi,
                                           G4double r0,
                                           G4double tanstereo,
                                           G4double zmin,
                                           G4double zmax);

  private:

    /**
     * Helper function for use by TriangulatePolygon().
     */
    static G4bool CheckSnip(const G4TwoVectorList& contour,
                            G4int a, G4int b, G4int c,
                            G4int n, const G4int* V);

    /**
     * Complete Elliptic Integral of the Second Kind.
     */
    static G4double comp_ellint_2(G4double e);
};

#endif
