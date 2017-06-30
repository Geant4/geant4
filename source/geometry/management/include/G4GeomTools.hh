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
//
// $Id: $
//
//
// --------------------------------------------------------------------
// class G4GeomTools
//
// Class description:
//
// A collection of utilities which can be helpfull for a wide range
// of geometry-related tasks

// History:
//
// 10.10.2016 E.Tcherniaev: Initial version.
// --------------------------------------------------------------------

#ifndef G4GEOMTOOLS_HH
#define G4GEOMTOOLS_HH

#include <vector>
#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"

typedef std::vector<G4TwoVector>   G4TwoVectorList;
typedef std::vector<G4ThreeVector> G4ThreeVectorList;

class G4GeomTools
{
  public:

  // ==================================================================
  //   2D Utilities 
  // ------------------------------------------------------------------

  static G4double TriangleArea(G4double Ax, G4double Ay,
                               G4double Bx, G4double By,
                               G4double Cx, G4double Cy);

  static G4double TriangleArea(const G4TwoVector& A,
                               const G4TwoVector& B,
                               const G4TwoVector& C);
    // Calcuate area of 2D triangle, return value is positive if
    // vertices of the triangle are given in anticlockwise order,
    // otherwise it is negative

  static G4double QuadArea(const G4TwoVector& A,
                           const G4TwoVector& B,
                           const G4TwoVector& C,
                           const G4TwoVector& D);
    // Calcuate area of 2D quadrilateral, return value is positive if
    // vertices of the quadrilateral are given in anticlockwise order,
    // otherwise it is negative

  static G4double PolygonArea(const G4TwoVectorList& polygon);
    // Calcuate area of 2D polygon, return value is positive if
    // vertices of the polygon are defined in anticlockwise order,
    // otherwise it is negative

  static G4bool PointInTriangle(G4double Px, G4double Py,
                                G4double Ax, G4double Ay,
                                G4double Bx, G4double By,
                                G4double Cx, G4double Cy);
    // Decide if point (Px,Py) is inside triangle (Ax,Ay)(Bx,By)(Cx,Cy)

  static G4bool PointInTriangle(const G4TwoVector& P,
                                const G4TwoVector& A,
                                const G4TwoVector& B,
                                const G4TwoVector& C);
    // Decide if point P is inside triangle ABC

  /*
  static G4bool PointInPolygon(const G4TwoVector& P
                               const G4TwoVectorList& Polygon);
    // Decide if point P is inside Polygon
  */

  static G4bool IsConvex(const G4TwoVectorList& polygon);
    // Decide if 2D polygon is convex, i.e. all internal angles are
    // less than pi

  static G4bool TriangulatePolygon(const G4TwoVectorList& polygon,
                                         G4TwoVectorList& result);

  static G4bool TriangulatePolygon(const G4TwoVectorList& polygon,
                                         std::vector<G4int>& result);
    // Simple implementation of "ear clipping" algorithm for
    // triangulation of a simple contour/polygon, it places results
    // in a std::vector as triplets of vertices. If triangulation
    // is sucsessfull then the function returns true, otherwise false

  static void RemoveRedundantVertices(G4TwoVectorList& polygon,
                                      std::vector<G4int>& iout,
                                      G4double tolerance = 0);
    // Remove collinear and coincident points from 2D polygon.
    // Indices of removed points are available in iout. 

  static G4bool DiskExtent(G4double rmin, G4double rmax,
                           G4double startPhi, G4double delPhi,
                           G4TwoVector& pmin, G4TwoVector& pmax); 
    // Calculate bounding rectangle of a disk sector,
    // it returns false if input parameters do not meet the following:
    //   rmin   >= 0
    //   rmax   >  rmin + kCarTolerance
    //   delPhi >  0 + kCarTolerance

  static void DiskExtent(G4double rmin, G4double rmax,
                         G4double sinPhiStart, G4double cosPhiStart,
                         G4double sinPhiEnd, G4double cosPhiEnd,
                         G4TwoVector& pmin, G4TwoVector& pmax); 
    // Calculate bounding rectangle of a disk sector,
    // faster version without check of parameters

  static G4double EllipsePerimeter(G4double a,
                                   G4double b);
    // Compute the circumference (perimeter) of an ellipse

  static G4double EllipticConeLateralArea(G4double a,
                                          G4double b,
                                          G4double h);
    // Compute the lateral surface area of an elliptic cone

  // ==================================================================
  //   3D Utilities 
  // ------------------------------------------------------------------

  static G4ThreeVector TriangleAreaNormal(const G4ThreeVector& A,
                                          const G4ThreeVector& B,
                                          const G4ThreeVector& C);
    // Find the normal to the plane of 3D triangle ABC,
    // length of the normal is equal to the area of the triangle

  static G4ThreeVector QuadAreaNormal(const G4ThreeVector& A,
                                      const G4ThreeVector& B,
                                      const G4ThreeVector& C,
                                      const G4ThreeVector& D);
    // Find normal to the plane of 3D quadrilateral ABCD,
    // length of the normal is equal to the area of the quadrilateral

  static G4ThreeVector PolygonAreaNormal(const G4ThreeVectorList& polygon);
    // Find normal to the plane of 3D polygon
    // length of the normal is equal to the area of the polygon

  /*
  static G4bool IsPlanar(const G4ThreeVector& A,
                         const G4ThreeVector& B,
                         const G4ThreeVector& C,
                         const G4ThreeVector& D);
    // Decide if 3D quadrilateral ABCD is planar

  static G4bool IsPlanar(const G4ThreeVectorList& polygon
                         const G4ThreeVector& normal);
    // Decide if 3D polygon is planar
 */

  static G4double DistancePointSegment(const G4ThreeVector& P,
                                       const G4ThreeVector& A,
                                       const G4ThreeVector& B);
    // Calculate distance between point P and line segment AB in 3D

  static G4ThreeVector ClosestPointOnSegment(const G4ThreeVector& P,
                                             const G4ThreeVector& A,
                                             const G4ThreeVector& B);
    // Find point on 3D line segment AB closest to point P

  static G4ThreeVector ClosestPointOnTriangle(const G4ThreeVector& P,
                                              const G4ThreeVector& A,
                                              const G4ThreeVector& B,
                                              const G4ThreeVector& C);
    // Find point on 3D triangle ABC closest to point P

  static G4bool SphereExtent(G4double rmin, G4double rmax,
                             G4double startTheta, G4double delTheta,
                             G4double startPhi, G4double delPhi,
                             G4ThreeVector& pmin, G4ThreeVector& pmax); 
    // Calculate bounding box of a spherical sector,
    // it returns false if input parameters do not meet the following:
    //   rmin       >= 0
    //   rmax       >  rmin + kCarTolerance
    //   startTheta >= 0 && <= pi;
    //   delTheta   >  0 + kCarTolerance
    //   delPhi     >  0 + kCarTolerance

  private:

  static G4bool CheckSnip(const G4TwoVectorList& contour,
                          G4int a, G4int b, G4int c,
                          G4int n, const G4int* V);
    // Helper function for TriangulatePolygon()

  static G4double comp_ellint_2(G4double e);
    // Complete Elliptic Integral of the Second Kind
};

#endif // G4GEOMTOOLS_HH
