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
// class G4GeomTools implementation
//
// Author: evgueni.tcherniaev@cern.ch
//
// 10.10.2016 E.Tcherniaev: initial version.
// --------------------------------------------------------------------

#include "G4GeomTools.hh"

#include "geomdefs.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryTolerance.hh"

///////////////////////////////////////////////////////////////////////
//
// Calculate area of a triangle in 2D

G4double G4GeomTools::TriangleArea(G4double Ax, G4double Ay,
                                   G4double Bx, G4double By,
                                   G4double Cx, G4double Cy)
{
  return ((Bx-Ax)*(Cy-Ay) - (By-Ay)*(Cx-Ax))*0.5;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate area of a triangle in 2D

G4double G4GeomTools::TriangleArea(const G4TwoVector& A,
                                   const G4TwoVector& B,
                                   const G4TwoVector& C)
{
  G4double Ax = A.x(), Ay = A.y(); 
  return ((B.x()-Ax)*(C.y()-Ay) - (B.y()-Ay)*(C.x()-Ax))*0.5;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate area of a quadrilateral in 2D

G4double G4GeomTools::QuadArea(const G4TwoVector& A,
                               const G4TwoVector& B,
                               const G4TwoVector& C,
                               const G4TwoVector& D)
{
  return ((C.x()-A.x())*(D.y()-B.y()) - (C.y()-A.y())*(D.x()-B.x()))*0.5;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate area of a polygon in 2D

G4double G4GeomTools::PolygonArea(const G4TwoVectorList& p)
{
  G4int n = p.size();
  if (n < 3) return 0; // degerate polygon
  G4double area = p[n-1].x()*p[0].y() - p[0].x()*p[n-1].y();
  for(G4int i=1; i<n; ++i)
  {
    area += p[i-1].x()*p[i].y() - p[i].x()*p[i-1].y();
  }
  return area*0.5;
}

///////////////////////////////////////////////////////////////////////
//
// Point inside 2D triangle

G4bool G4GeomTools::PointInTriangle(G4double Ax, G4double Ay,
                                    G4double Bx, G4double By,
                                    G4double Cx, G4double Cy,
                                    G4double Px, G4double Py)

{
  if ((Bx-Ax)*(Cy-Ay) - (By-Ay)*(Cx-Ax) > 0.)
  {
    if ((Ax-Cx)*(Py-Cy) - (Ay-Cy)*(Px-Cx) < 0.) return false;
    if ((Bx-Ax)*(Py-Ay) - (By-Ay)*(Px-Ax) < 0.) return false;
    if ((Cx-Bx)*(Py-By) - (Cy-By)*(Px-Bx) < 0.) return false;
  }
  else
  {
    if ((Ax-Cx)*(Py-Cy) - (Ay-Cy)*(Px-Cx) > 0.) return false;
    if ((Bx-Ax)*(Py-Ay) - (By-Ay)*(Px-Ax) > 0.) return false;
    if ((Cx-Bx)*(Py-By) - (Cy-By)*(Px-Bx) > 0.) return false;
  }
  return true;
}

///////////////////////////////////////////////////////////////////////
//
// Point inside 2D triangle

G4bool G4GeomTools::PointInTriangle(const G4TwoVector& A,
                                    const G4TwoVector& B,
                                    const G4TwoVector& C,
                                    const G4TwoVector& P)
{
  G4double Ax = A.x(), Ay = A.y();
  G4double Bx = B.x(), By = B.y();
  G4double Cx = C.x(), Cy = C.y();
  G4double Px = P.x(), Py = P.y();
  if ((Bx-Ax)*(Cy-Ay) - (By-Ay)*(Cx-Ax) > 0.)
  {
    if ((Ax-Cx)*(Py-Cy) - (Ay-Cy)*(Px-Cx) < 0.) return false;
    if ((Bx-Ax)*(Py-Ay) - (By-Ay)*(Px-Ax) < 0.) return false;
    if ((Cx-Bx)*(Py-By) - (Cy-By)*(Px-Bx) < 0.) return false;
  }
  else
  {
    if ((Ax-Cx)*(Py-Cy) - (Ay-Cy)*(Px-Cx) > 0.) return false;
    if ((Bx-Ax)*(Py-Ay) - (By-Ay)*(Px-Ax) > 0.) return false;
    if ((Cx-Bx)*(Py-By) - (Cy-By)*(Px-Bx) > 0.) return false;
  }
  return true;
}

///////////////////////////////////////////////////////////////////////
//
// Detemine whether 2D polygon is convex or not

G4bool G4GeomTools::IsConvex(const G4TwoVectorList& polygon)
{
  static const G4double kCarTolerance =
         G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  G4bool gotNegative = false;
  G4bool gotPositive = false;
  G4int n = polygon.size();
  if (n <= 0) return false;
  for (G4int icur=0; icur<n; ++icur)
  {
    G4int iprev = (icur ==   0) ? n-1 : icur-1;
    G4int inext = (icur == n-1) ?   0 : icur+1;
    G4TwoVector e1 = polygon[icur]  - polygon[iprev];
    G4TwoVector e2 = polygon[inext] - polygon[icur];
    G4double cross = e1.x()*e2.y() - e1.y()*e2.x();
    if (std::abs(cross) < kCarTolerance) return false;
    if (cross <  0) gotNegative = true;
    if (cross >  0) gotPositive = true;
    if (gotNegative && gotPositive) return false;
  }
  return true;
}

///////////////////////////////////////////////////////////////////////
//
// Triangulate simple polygon

G4bool G4GeomTools::TriangulatePolygon(const G4TwoVectorList& polygon,
                                             G4TwoVectorList& result)
{
  result.resize(0);
  std::vector<G4int> triangles;
  G4bool reply = TriangulatePolygon(polygon,triangles);

  G4int n = triangles.size();
  for (G4int i=0; i<n; ++i) result.push_back(polygon[triangles[i]]);
  return reply;
}

///////////////////////////////////////////////////////////////////////
//
// Triangulation of a simple polygon by "ear clipping"

G4bool G4GeomTools::TriangulatePolygon(const G4TwoVectorList& polygon,
                                             std::vector<G4int>& result)
{
  result.resize(0);

  // allocate and initialize list of Vertices in polygon
  //
  G4int n = polygon.size();
  if (n < 3) return false;

  // we want a counter-clockwise polygon in V
  // 
  G4double area = G4GeomTools::PolygonArea(polygon);
  G4int* V = new G4int[n];
  if (area > 0.)
    for (G4int i=0; i<n; ++i) V[i] = i;
  else
    for (G4int i=0; i<n; ++i) V[i] = (n-1)-i;

  //  Triangulation: remove nv-2 Vertices, creating 1 triangle every time
  // 
  G4int nv = n;
  G4int count = 2*nv; // error detection counter
  for(G4int b=nv-1; nv>2; )
  {
    // ERROR: if we loop, it is probably a non-simple polygon
    if ((count--) <= 0)
    {
      delete[] V;
      if (area < 0.) std::reverse(result.begin(),result.end());
      return false; 
    }

    // three consecutive vertices in current polygon, <a,b,c>
    G4int a = (b   < nv) ? b   : 0; // previous
          b = (a+1 < nv) ? a+1 : 0; // current
    G4int c = (b+1 < nv) ? b+1 : 0; // next

    if (CheckSnip(polygon, a,b,c, nv,V))
    {
      // output Triangle
      result.push_back(V[a]);
      result.push_back(V[b]);
      result.push_back(V[c]);

      // remove vertex b from remaining polygon
      nv--;
      for(G4int i=b; i<nv; ++i) V[i] = V[i+1];

      count = 2*nv; // resest error detection counter
    }
  }
  delete[] V;
  if (area < 0.) std::reverse(result.begin(),result.end());
  return true;
}

///////////////////////////////////////////////////////////////////////
//
// Helper function for "ear clipping" polygon triangulation.
// Check for a valid snip

G4bool G4GeomTools::CheckSnip(const G4TwoVectorList& contour,
                              G4int a, G4int b, G4int c,
                              G4int n, const G4int* V)
{
  static const G4double kCarTolerance =
         G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  // check orientation of Triangle
  G4double Ax = contour[V[a]].x(), Ay = contour[V[a]].y();
  G4double Bx = contour[V[b]].x(), By = contour[V[b]].y();
  G4double Cx = contour[V[c]].x(), Cy = contour[V[c]].y();
  if ((Bx-Ax)*(Cy-Ay) - (By-Ay)*(Cx-Ax) < kCarTolerance) return false;
  
  // check that there is no point inside Triangle
  G4double xmin = std::min(std::min(Ax,Bx),Cx);
  G4double xmax = std::max(std::max(Ax,Bx),Cx);
  G4double ymin = std::min(std::min(Ay,By),Cy);
  G4double ymax = std::max(std::max(Ay,By),Cy);
  for (G4int i=0; i<n; ++i)
  {
    if((i == a) || (i == b) || (i == c)) continue;
    G4double Px = contour[V[i]].x();
    if (Px < xmin || Px > xmax) continue;
    G4double Py = contour[V[i]].y();
    if (Py < ymin || Py > ymax) continue;
    if (PointInTriangle(Ax,Ay,Bx,By,Cx,Cy,Px,Py)) return false;
  }
  return true;
}

///////////////////////////////////////////////////////////////////////
//
// Remove collinear and coincident points from 2D polygon

void G4GeomTools::RemoveRedundantVertices(G4TwoVectorList& polygon,
                                          std::vector<G4int>& iout,
                                          G4double tolerance) 
{
  iout.resize(0);
  // set tolerance squared
  G4double delta = tolerance*tolerance;
  // set special value to mark vertices for removal
  G4double removeIt = kInfinity;

  G4int nv = polygon.size();

  // Main loop: check every three consecutive points, if the points
  // are collinear then mark middle point for removal
  //
  G4int icur = 0, iprev = 0, inext = 0, nout = 0; 
  for (G4int i=0; i<nv; ++i)
  {
    icur = i;                    // index of current point

    for (G4int k=1; k<nv+1; ++k) // set index of previous point
    {
      iprev = icur - k;
      if (iprev < 0) iprev += nv;
      if (polygon[iprev].x() != removeIt) break;
    }

    for (G4int k=1; k<nv+1; ++k) // set index of next point
    {
      inext = icur + k;
      if (inext >= nv) inext -= nv;
      if (polygon[inext].x() != removeIt) break;
    }

    if (iprev == inext) break;   // degenerate polygon, stop

    // Calculate parameters of triangle (iprev->icur->inext),
    // if triangle is too small or too narrow then mark current
    // point for removal
    G4TwoVector e1 = polygon[iprev] - polygon[icur];
    G4TwoVector e2 = polygon[inext] - polygon[icur];

    // Check length of edges, then check height of the triangle 
    G4double leng1 = e1.mag2();
    G4double leng2 = e2.mag2();
    G4double leng3 = (e2-e1).mag2();
    if (leng1 <= delta || leng2 <= delta || leng3 <= delta)
    {
      polygon[icur].setX(removeIt); nout++;
    }
    else
    {
      G4double lmax = std::max(std::max(leng1,leng2),leng3);
      G4double area = std::abs(e1.x()*e2.y()-e1.y()*e2.x())*0.5;
      if (area/std::sqrt(lmax) <= std::abs(tolerance))
      {
        polygon[icur].setX(removeIt); nout++;
      }
    }
  }

  // Remove marked points
  //
  icur = 0;
  if (nv - nout < 3)           // degenerate polygon, remove all points
  {
    for (G4int i=0; i<nv; ++i) iout.push_back(i);
    polygon.resize(0);
    nv = 0;
  }
  for (G4int i=0; i<nv; ++i) // move points, if required
  {
    if (polygon[i].x() != removeIt)
      polygon[icur++] = polygon[i];
    else
      iout.push_back(i);
  }
  if (icur < nv) polygon.resize(icur);
  return;
}

///////////////////////////////////////////////////////////////////////
//
// Find bounding rectangle of a disk sector

G4bool G4GeomTools::DiskExtent(G4double rmin, G4double rmax,
                               G4double startPhi, G4double delPhi,
                               G4TwoVector& pmin, G4TwoVector& pmax)
{
  static const G4double kCarTolerance =
         G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  // check parameters
  //
  pmin.set(0,0);
  pmax.set(0,0);
  if (rmin   <  0)                    return false;
  if (rmax   <= rmin + kCarTolerance) return false;
  if (delPhi <= 0    + kCarTolerance) return false;

  // calculate extent
  //
  pmin.set(-rmax,-rmax);
  pmax.set( rmax, rmax);
  if (delPhi >= CLHEP::twopi) return true;

  DiskExtent(rmin,rmax,
             std::sin(startPhi),std::cos(startPhi),
             std::sin(startPhi+delPhi),std::cos(startPhi+delPhi),
             pmin,pmax);
  return true;
}

///////////////////////////////////////////////////////////////////////
//
// Find bounding rectangle of a disk sector, fast version.
// No check of parameters !!!

void G4GeomTools::DiskExtent(G4double rmin, G4double rmax,
                             G4double sinStart, G4double cosStart,
                             G4double sinEnd, G4double cosEnd,
                             G4TwoVector& pmin, G4TwoVector& pmax)
{
  static const G4double kCarTolerance =
         G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  // check if 360 degrees
  //
  pmin.set(-rmax,-rmax);
  pmax.set( rmax, rmax);

  if (std::abs(sinEnd-sinStart) < kCarTolerance && 
      std::abs(cosEnd-cosStart) < kCarTolerance) return;

  // get start and end quadrants
  //
  //      1 | 0
  //     ---+--- 
  //      3 | 2
  //
  G4int icase = (cosEnd < 0) ? 1 : 0;
  if (sinEnd   < 0) icase += 2;
  if (cosStart < 0) icase += 4;
  if (sinStart < 0) icase += 8;

  switch (icase)
  {
  // start quadrant 0
  case  0:                                 // start->end : 0->0
    if (sinEnd < sinStart) break;
    pmin.set(rmin*cosEnd,rmin*sinStart);
    pmax.set(rmax*cosStart,rmax*sinEnd  );
    break;
  case  1:                                 // start->end : 0->1
    pmin.set(rmax*cosEnd,std::min(rmin*sinStart,rmin*sinEnd));
    pmax.set(rmax*cosStart,rmax  );
    break;
  case  2:                                 // start->end : 0->2
    pmin.set(-rmax,-rmax);
    pmax.set(std::max(rmax*cosStart,rmax*cosEnd),rmax);
    break;
  case  3:                                 // start->end : 0->3
    pmin.set(-rmax,rmax*sinEnd);
    pmax.set(rmax*cosStart,rmax);
    break;
  // start quadrant 1
  case  4:                                 // start->end : 1->0
    pmin.set(-rmax,-rmax);
    pmax.set(rmax,std::max(rmax*sinStart,rmax*sinEnd));
    break;
  case  5:                                 // start->end : 1->1
    if (sinEnd > sinStart) break;
    pmin.set(rmax*cosEnd,rmin*sinEnd  );
    pmax.set(rmin*cosStart,rmax*sinStart);
    break;
  case  6:                                 // start->end : 1->2
    pmin.set(-rmax,-rmax);
    pmax.set(rmax*cosEnd,rmax*sinStart);
    break;
  case  7:                                 // start->end : 1->3
    pmin.set(-rmax,rmax*sinEnd);
    pmax.set(std::max(rmin*cosStart,rmin*cosEnd),rmax*sinStart);
    break;
  // start quadrant 2
  case  8:                                 // start->end : 2->0
    pmin.set(std::min(rmin*cosStart,rmin*cosEnd),rmax*sinStart);
    pmax.set(rmax,rmax*sinEnd);
    break;
  case  9:                                 // start->end : 2->1
    pmin.set(rmax*cosEnd,rmax*sinStart);
    pmax.set(rmax,rmax);
    break;
  case 10:                                 // start->end : 2->2
    if (sinEnd < sinStart) break;
    pmin.set(rmin*cosStart,rmax*sinStart);
    pmax.set(rmax*cosEnd,rmin*sinEnd  );
    break;
  case 11:                                 // start->end : 2->3
    pmin.set(-rmax,std::min(rmax*sinStart,rmax*sinEnd));
    pmax.set(rmax,rmax);
    break;
  // start quadrant 3
  case 12:                                 // start->end : 3->0
    pmin.set(rmax*cosStart,-rmax);
    pmax.set(rmax,rmax*sinEnd);
    break;
  case 13:                                 // start->end : 3->1
    pmin.set(std::min(rmax*cosStart,rmax*cosEnd),-rmax);
    pmax.set(rmax,rmax);
    break;
  case 14:                                 // start->end : 3->2
    pmin.set(rmax*cosStart,-rmax);
    pmax.set(rmax*cosEnd,std::max(rmin*sinStart,rmin*sinEnd));
    break;
  case 15:                                 // start->end : 3->3
    if (sinEnd > sinStart) break;
    pmin.set(rmax*cosStart,rmax*sinEnd);
    pmax.set(rmin*cosEnd,rmin*sinStart);
    break;
  }
  return;
}

///////////////////////////////////////////////////////////////////////
//
// Compute the circumference (perimeter) of an ellipse

G4double G4GeomTools::EllipsePerimeter(G4double pA, G4double pB)
{
  G4double x = std::abs(pA);
  G4double y = std::abs(pB);
  G4double a = std::max(x,y);
  G4double b = std::min(x,y);
  G4double e = std::sqrt((1. - b/a)*(1. + b/a));
  return 4. * a * comp_ellint_2(e);
}

///////////////////////////////////////////////////////////////////////
//
// Compute the lateral surface area of an elliptic cone

G4double G4GeomTools::EllipticConeLateralArea(G4double pA,
                                              G4double pB,
                                              G4double pH)
{
  G4double x = std::abs(pA);
  G4double y = std::abs(pB);
  G4double h = std::abs(pH);
  G4double a = std::max(x,y);
  G4double b = std::min(x,y);
  G4double e = std::sqrt((1. - b/a)*(1. + b/a)) / std::hypot(1.,b/h);
  return 2. * a * std::hypot(b,h) * comp_ellint_2(e);
}

///////////////////////////////////////////////////////////////////////
//
// Compute Elliptical Integral of the Second Kind
//
// The algorithm is based upon Carlson B.C., "Computation of real
// or complex elliptic integrals", Numerical Algorithms,
// Volume 10, Issue 1, 1995 (see equations 2.36 - 2.39)
//
// The code was adopted from C code at:
// http://paulbourke.net/geometry/ellipsecirc/

G4double G4GeomTools::comp_ellint_2(G4double e)
{
  const G4double eps = 1. / 134217728.; // 1 / 2^27

  G4double a = 1.;
  G4double b = std::sqrt((1. - e)*(1. + e));
  if (b == 1.) return CLHEP::halfpi;
  if (b == 0.) return 1.;

  G4double x = 1.;
  G4double y = b;
  G4double S = 0.;
  G4double M = 1.;
  while (x - y > eps*y) {
    G4double tmp = (x + y) * 0.5;
    y = std::sqrt(x*y);
    x = tmp;
    M += M;
    S += M * (x - y)*(x - y);
  }
  return 0.5 * CLHEP::halfpi * ((a + b)*(a + b) - S) / (x + y);
}

///////////////////////////////////////////////////////////////////////
//
// Calcuate area of a triangle in 3D

G4ThreeVector G4GeomTools::TriangleAreaNormal(const G4ThreeVector& A,
                                              const G4ThreeVector& B,
                                              const G4ThreeVector& C)
{
  return ((B-A).cross(C-A))*0.5;
}

///////////////////////////////////////////////////////////////////////
//
// Calcuate area of a quadrilateral in 3D

G4ThreeVector G4GeomTools::QuadAreaNormal(const G4ThreeVector& A,
                                          const G4ThreeVector& B,
                                          const G4ThreeVector& C,
                                          const G4ThreeVector& D)
{
  return ((C-A).cross(D-B))*0.5;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate area of a polygon in 3D

G4ThreeVector G4GeomTools::PolygonAreaNormal(const G4ThreeVectorList& p)
{
  G4int n = p.size();
  if (n < 3) return G4ThreeVector(0,0,0); // degerate polygon
  G4ThreeVector normal = p[n-1].cross(p[0]);
  for(G4int i=1; i<n; ++i)
  {
    normal += p[i-1].cross(p[i]);
  }
  return normal*0.5;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate distance between point P and line segment AB in 3D

G4double G4GeomTools::DistancePointSegment(const G4ThreeVector& P,
                                           const G4ThreeVector& A,
                                           const G4ThreeVector& B)
{
  G4ThreeVector AP = P - A;
  G4ThreeVector AB = B - A;

  G4double u = AP.dot(AB);
  if (u <= 0) return AP.mag();       // closest point is A

  G4double len2 = AB.mag2();
  if (u >= len2) return (B-P).mag(); // closest point is B

  return ((u/len2)*AB - AP).mag();   // distance to line
}

///////////////////////////////////////////////////////////////////////
//
// Find closest point on line segment in 3D

G4ThreeVector
G4GeomTools::ClosestPointOnSegment(const G4ThreeVector& P,
                                   const G4ThreeVector& A,
                                   const G4ThreeVector& B)
{
  G4ThreeVector AP = P - A;
  G4ThreeVector AB = B - A;

  G4double u = AP.dot(AB);
  if (u <= 0) return A;      // closest point is A

  G4double len2 = AB.mag2();
  if (u >= len2) return B;   // closest point is B

  G4double t = u/len2;
  return A + t*AB;           // closest point on segment
}

///////////////////////////////////////////////////////////////////////
//
// Find closest point on triangle in 3D.
//
// The implementation is based on the algorithm published in
// "Geometric Tools for Computer Graphics", Philip J Scheider and
// David H Eberly, Elsevier Science (USA), 2003.
// 
// The algorithm is also available at:
// http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf

G4ThreeVector
G4GeomTools::ClosestPointOnTriangle(const G4ThreeVector& P,
                                    const G4ThreeVector& A,
                                    const G4ThreeVector& B,
                                    const G4ThreeVector& C)
{
  G4ThreeVector diff  = A - P;
  G4ThreeVector edge0 = B - A;
  G4ThreeVector edge1 = C - A;

  G4double a = edge0.mag2();
  G4double b = edge0.dot(edge1);
  G4double c = edge1.mag2();
  G4double d = diff.dot(edge0);
  G4double e = diff.dot(edge1);

  G4double det = a*c - b*b;
  G4double t0  = b*e - c*d;
  G4double t1  = b*d - a*e;

  /*
             ^ t1
         \ 2 |
          \  |
           \ |     regions
            \|
             C
             |\
         3   | \   1
             |  \
             | 0 \
             |    \
        ---- A --- B ----> t0
             |      \
         4   |   5   \   6
             |        \
  */

  G4int region = -1;
  if (t0+t1 <= det)
    region = (t0 < 0) ? ((t1 < 0) ? 4 : 3) : ((t1 < 0) ? 5 : 0);
  else
    region = (t0 < 0) ? 2 : ((t1 < 0) ? 6 : 1);

  switch (region)
  {
  case 0:  // interior of triangle
    {
      G4double invDet = 1./det;
      return A + (t0*invDet)*edge0 + (t1*invDet)*edge1;
    }
  case 1:  // edge BC
    {
      G4double numer = c + e - b - d; 
      if (numer <= 0) return C;
      G4double denom = a - 2*b + c;
      return (numer >= denom) ? B : C + (numer/denom)*(edge0-edge1);
    }
  case 2:  // edge AC or BC
    {
      G4double tmp0 = b + d;
      G4double tmp1 = c + e;
      if (tmp1 > tmp0)
      {
        G4double numer = tmp1 - tmp0;
        G4double denom = a - 2*b + c;
        return (numer >= denom) ? B : C + (numer/denom)*(edge0-edge1);
      }
      // same:  (e >= 0) ? A : ((-e >= c) ? C : A + (-e/c)*edge1)
      return (tmp1 <= 0) ? C : (( e >= 0) ? A : A + (-e/c)*edge1);
    }
  case 3:  // edge AC
    return (e >= 0) ? A : ((-e >= c) ? C : A + (-e/c)*edge1);

  case 4:  // edge AB or AC
    if (d < 0)      return (-d >= a) ? B : A + (-d/a)*edge0;
    return (e >= 0) ? A : ((-e >= c) ? C : A + (-e/c)*edge1);

  case 5:  // edge AB
    return (d >= 0) ? A : ((-d >= a) ? B : A + (-d/a)*edge0);

  case 6:  // edge AB or BC
    {
      G4double tmp0 = b + e;
      G4double tmp1 = a + d;
      if (tmp1 > tmp0)
      {
        G4double numer = tmp1 - tmp0;
        G4double denom = a - 2*b + c;
        return (numer >= denom) ? C : B + (numer/denom)*(edge1-edge0);
      }
      // same:  (d >= 0) ? A : ((-d >= a) ? B : A + (-d/a)*edge0)
      return (tmp1 <= 0) ? B : (( d >= 0) ? A : A + (-d/a)*edge0);
    }
  default: // impossible case 
    return G4ThreeVector(kInfinity,kInfinity,kInfinity);
  }
}

///////////////////////////////////////////////////////////////////////
//
// Calculate bounding box of a spherical sector

G4bool
G4GeomTools::SphereExtent(G4double rmin, G4double rmax,
                          G4double startTheta, G4double delTheta,
                          G4double startPhi, G4double delPhi,
                          G4ThreeVector& pmin, G4ThreeVector& pmax)
{
  static const G4double kCarTolerance =
         G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  // check parameters
  //
  pmin.set(0,0,0);
  pmax.set(0,0,0);
  if (rmin     <  0)                    return false;
  if (rmax     <= rmin + kCarTolerance) return false;
  if (delTheta <= 0    + kCarTolerance) return false;
  if (delPhi   <= 0    + kCarTolerance) return false;

  G4double stheta = startTheta;
  G4double dtheta = delTheta;
  if (stheta < 0 && stheta > CLHEP::pi) return false;
  if (stheta + dtheta > CLHEP::pi)      dtheta = CLHEP::pi - stheta;
  if (dtheta <= 0 + kCarTolerance)      return false;

  // calculate extent
  //
  pmin.set(-rmax,-rmax,-rmax);
  pmax.set( rmax, rmax, rmax);
  if (dtheta >= CLHEP::pi && delPhi >= CLHEP::twopi) return true;

  G4double etheta   = stheta + dtheta;
  G4double sinStart = std::sin(stheta);
  G4double cosStart = std::cos(stheta);
  G4double sinEnd   = std::sin(etheta);
  G4double cosEnd   = std::cos(etheta);

  G4double rhomin = rmin*std::min(sinStart,sinEnd);
  G4double rhomax = rmax;
  if (stheta > CLHEP::halfpi) rhomax = rmax*sinStart;
  if (etheta < CLHEP::halfpi) rhomax = rmax*sinEnd;

  G4TwoVector xymin,xymax;
  DiskExtent(rhomin,rhomax,
             std::sin(startPhi),std::cos(startPhi),
             std::sin(startPhi+delPhi),std::cos(startPhi+delPhi),
             xymin,xymax);

  G4double zmin = std::min(rmin*cosEnd,rmax*cosEnd);
  G4double zmax = std::max(rmin*cosStart,rmax*cosStart);
  pmin.set(xymin.x(),xymin.y(),zmin);
  pmax.set(xymax.x(),xymax.y(),zmax);
  return true;
}
