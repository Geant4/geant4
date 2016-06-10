//
// ********************************************************************
// * This Software is part of the AIDA Unified Solids Library package *
// * See: https://aidasoft.web.cern.ch/USolids                        *
// ********************************************************************
//
// $Id:$
//
// --------------------------------------------------------------------
//
// UPolyPhiFace
//
// 19.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include <string>
#include <cmath>
#include <sstream>
#include "UPolyPhiFace.hh"
#include "UReduciblePolygon.hh"
#include "UVector2.hh"

//
// Constructor
//
// Points r,z should be supplied in clockwise order in r,z. For example:
//
//                [1]---------[2]        ^ R
//                 |           |          |
//                 |           |          +--> z
//                [0]---------[3]
//
UPolyPhiFace::UPolyPhiFace(const UReduciblePolygon* rz,
                           double phi,
                           double deltaPhi,
                           double phiOther)
  : fSurfaceArea(0.), triangles(0)
{
  fTolerance = VUSolid::Tolerance();

  numEdges = rz->NumVertices();

  rMin = rz->Amin();
  rMax = rz->Amax();
  zMin = rz->Bmin();
  zMax = rz->Bmax();

  //
  // Is this the "starting" phi edge of the two?
  //
  bool start = (phiOther > phi);

  //
  // Build radial vector
  //
  radial = UVector3(std::cos(phi), std::sin(phi), 0.0);

  //
  // Build normal
  //
  double zSign = start ? 1 : -1;
  normal = UVector3(zSign * radial.y(), -zSign * radial.x(), 0);

  //
  // Is allBehind?
  //
  allBehind = (zSign * (std::cos(phiOther) * radial.y() - std::sin(phiOther) * radial.x()) < 0);

  //
  // Adjacent edges
  //
  double midPhi = phi + (start ? +0.5 : -0.5) * deltaPhi;
  double cosMid = std::cos(midPhi),
         sinMid = std::sin(midPhi);

  //
  // Allocate corners
  //
  corners = new UPolyPhiFaceVertex[numEdges];
  //
  // Fill them
  //
  UReduciblePolygonIterator iterRZ(rz);

  UPolyPhiFaceVertex* corn = corners;
  UPolyPhiFaceVertex* helper = corners;

  iterRZ.Begin();
  do
  {
    corn->r = iterRZ.GetA();
    corn->z = iterRZ.GetB();
    corn->x = corn->r * radial.x();
    corn->y = corn->r * radial.y();

    // Add pointer on prev corner
    //
    if (corn == corners)
    {
      corn->prev = corners + numEdges - 1;
    }
    else
    {
      corn->prev = helper;
    }

    // Add pointer on next corner
    //
    if (corn < corners + numEdges - 1)
    {
      corn->next = corn + 1;
    }
    else
    {
      corn->next = corners;
    }

    helper = corn;
  }
  while (++corn, iterRZ.Next());

  //
  // Allocate edges
  //
  edges = new UPolyPhiFaceEdge[numEdges];

  //
  // Fill them
  //
  double rFact = std::cos(0.5 * deltaPhi);
  double rFactNormalize = 1.0 / std::sqrt(1.0 + rFact * rFact);

  UPolyPhiFaceVertex* prev = corners + numEdges - 1,
                      *here = corners;
  UPolyPhiFaceEdge*   edge = edges;
  do
  {
    UVector3 sideNorm;

    edge->v0 = prev;
    edge->v1 = here;

    double dr = here->r - prev->r,
           dz = here->z - prev->z;

    edge->length = std::sqrt(dr * dr + dz * dz);

    edge->tr = dr / edge->length;
    edge->tz = dz / edge->length;

    if ((here->r < DBL_MIN) && (prev->r < DBL_MIN))
    {
      //
      // Sigh! Always exceptions!
      // This edge runs at r==0, so its adjoing surface is not a
      // PolyconeSide or PolyhedraSide, but the opposite PolyPhiFace.
      //
      double zSignOther = start ? -1 : 1;
      sideNorm = UVector3(zSignOther * std::sin(phiOther),
                          -zSignOther * std::cos(phiOther), 0);
    }
    else
    {
      sideNorm = UVector3(edge->tz * cosMid,
                          edge->tz * sinMid,
                          -edge->tr * rFact);
      sideNorm *= rFactNormalize;
    }
    sideNorm += normal;

    edge->norm3D = sideNorm.Unit();
  }
  while (edge++, prev = here, ++here < corners + numEdges);

  //
  // Go back and fill in corner "normals"
  //
  UPolyPhiFaceEdge* prevEdge = edges + numEdges - 1;
  edge = edges;
  do
  {
    //
    // Calculate vertex 2D normals (on the phi surface)
    //
    double rPart = prevEdge->tr + edge->tr;
    double zPart = prevEdge->tz + edge->tz;
    double norm = std::sqrt(rPart * rPart + zPart * zPart);
    double rNorm = +zPart / norm;
    double zNorm = -rPart / norm;

    edge->v0->rNorm = rNorm;
    edge->v0->zNorm = zNorm;

    //
    // Calculate the 3D normals.
    //
    // Find the vector perpendicular to the z axis
    // that defines the plane that contains the vertex normal
    //
    UVector3 xyVector;

    if (edge->v0->r < DBL_MIN)
    {
      //
      // This is a vertex at r==0, which is a special
      // case. The normal we will construct lays in the
      // plane at the center of the phi opening.
      //
      // We also know that rNorm < 0
      //
      double zSignOther = start ? -1 : 1;
      UVector3 normalOther(zSignOther * std::sin(phiOther),
                           -zSignOther * std::cos(phiOther), 0);

      xyVector = - normal - normalOther;
    }
    else
    {
      //
      // This is a vertex at r > 0. The plane
      // is the average of the normal and the
      // normal of the adjacent phi face
      //
      xyVector = UVector3(cosMid, sinMid, 0);
      if (rNorm < 0)
        xyVector -= normal;
      else
        xyVector += normal;
    }

    //
    // Combine it with the r/z direction from the face
    //
    edge->v0->norm3D = rNorm * xyVector.Unit() + UVector3(0, 0, zNorm);
  }
  while (prevEdge = edge, ++edge < edges + numEdges);

  //
  // Build point on surface
  //
  double rAve = 0.5 * (rMax - rMin),
         zAve = 0.5 * (zMax - zMin);
  surface = UVector3(rAve * radial.x(), rAve * radial.y(), zAve);
}


//
// Diagnose
//
// Throw an exception if something is found inconsistent with
// the solid.
//
// For debugging purposes only
//
void UPolyPhiFace::Diagnose(VUSolid* owner)
{
  UPolyPhiFaceVertex*   corner = corners;
  do
  {
    UVector3 test(corner->x, corner->y, corner->z);
    test -= 1E-6 * corner->norm3D;

    if (owner->Inside(test) != VUSolid::eInside)
    {
      UUtils::Exception("UPolyPhiFace::Diagnose()", "GeomSolids0002",
                        UFatalError, 1, "Bad vertex normal found.");
    }
  }
  while (++corner < corners + numEdges);
}


//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
UPolyPhiFace::UPolyPhiFace(__void__&)
  : numEdges(0), edges(0), corners(0), rMin(0.), rMax(0.), zMin(0.), zMax(0.),
    allBehind(false), fTolerance(0.), fSurfaceArea(0.), triangles(0)
{
}


//
// Destructor
//
UPolyPhiFace::~UPolyPhiFace()
{
  delete [] edges;
  delete [] corners;
}


//
// Copy constructor
//
UPolyPhiFace::UPolyPhiFace(const UPolyPhiFace& source)
  : UVCSGface()
{
  CopyStuff(source);
}


//
// Assignment operator
//
UPolyPhiFace& UPolyPhiFace::operator=(const UPolyPhiFace& source)
{
  if (this == &source)
  {
    return *this;
  }

  delete [] edges;
  delete [] corners;

  CopyStuff(source);

  return *this;
}


//
// CopyStuff (protected)
//
void UPolyPhiFace::CopyStuff(const UPolyPhiFace& source)
{
  //
  // The simple stuff
  //
  numEdges  = source.numEdges;
  normal    = source.normal;
  radial    = source.radial;
  surface  = source.surface;
  rMin    = source.rMin;
  rMax    = source.rMax;
  zMin    = source.zMin;
  zMax    = source.zMax;
  allBehind = source.allBehind;
  triangles = 0;

  fTolerance = source.fTolerance;
  fSurfaceArea = source.fSurfaceArea;

  //
  // Corner dynamic array
  //
  corners = new UPolyPhiFaceVertex[numEdges];
  UPolyPhiFaceVertex* corn = corners,
                      *sourceCorn = source.corners;
  do
  {
    *corn = *sourceCorn;
  }
  while (++sourceCorn, ++corn < corners + numEdges);

  //
  // Edge dynamic array
  //
  edges = new UPolyPhiFaceEdge[numEdges];

  UPolyPhiFaceVertex* prev = corners + numEdges - 1,
                      *here = corners;
  UPolyPhiFaceEdge*   edge = edges,
                      *sourceEdge = source.edges;
  do
  {
    *edge = *sourceEdge;
    edge->v0 = prev;
    edge->v1 = here;
  }
  while (++sourceEdge, ++edge, prev = here, ++here < corners + numEdges);
}


//
// Intersect
//
bool UPolyPhiFace::Distance(const UVector3& p,
                            const UVector3& v,
                            bool outgoing,
                            double surfTolerance,
                            double& distance,
                            double& distFromSurface,
                            UVector3& aNormal,
                            bool& isAllBehind)
{
  double normSign = outgoing ? +1 : -1;

  //
  // These don't change
  //
  isAllBehind = allBehind;
  aNormal = normal;

  //
  // Correct normal? Here we have straight sides, and can safely ignore
  // intersections where the Dot product with the normal is zero.
  //
  double dotProd = normSign * normal.Dot(v);

  if (dotProd <= 0) return false;

  //
  // Calculate distance to surface. If the side is too far
  // behind the point, we must reject it.
  //
  UVector3 ps = p - surface;
  distFromSurface = -normSign * ps.Dot(normal);

  if (distFromSurface < -surfTolerance) return false;

  //
  // Calculate precise distance to intersection with the side
  // (along the trajectory, not normal to the surface)
  //
  distance = distFromSurface / dotProd;

  //
  // Calculate intersection point in r,z
  //
  UVector3 ip = p + distance * v;

  double r = radial.Dot(ip);

  //
  // And is it inside the r/z extent?
  //
  return InsideEdgesExact(r, ip.z(), normSign, p, v);
}


//
// Distance
//
double UPolyPhiFace::Safety(const UVector3& p, bool outgoing)
{
  double normSign = outgoing ? +1 : -1;
  //
  // Correct normal?
  //
  UVector3 ps = p - surface;
  double distPhi = -normSign * normal.Dot(ps);

  if (distPhi < -0.5 * VUSolid::Tolerance())
    return UUtils::kInfinity;
  else if (distPhi < 0)
    distPhi = 0.0;

  //
  // Calculate projected point in r,z
  //
  double r = radial.Dot(p);

  //
  // Are we inside the face?
  //
  double distRZ2;

  if (InsideEdges(r, p.z(), &distRZ2, 0))
  {
    //
    // Yup, answer is just distPhi
    //
    return distPhi;
  }
  else
  {
    //
    // Nope. Penalize by distance out
    //
    return std::sqrt(distPhi * distPhi + distRZ2);
  }
}


//
// Inside
//
VUSolid::EnumInside UPolyPhiFace::Inside(const UVector3& p,
                                         double tolerance,
                                         double* bestDistance)
{
  //
  // Get distance along phi, which if negative means the point
  // is nominally inside the shape.
  //
  UVector3 ps = p - surface;
  double distPhi = normal.Dot(ps);

  //
  // Calculate projected point in r,z
  //
  double r = radial.Dot(p);

  //
  // Are we inside the face?
  //
  double distRZ2;
  UPolyPhiFaceVertex* base3Dnorm;
  UVector3*      head3Dnorm;

  if (InsideEdges(r, p.z(), &distRZ2, &base3Dnorm, &head3Dnorm))
  {
    //
    // Looks like we're inside. Distance is distance in phi.
    //
    *bestDistance = std::fabs(distPhi);

    //
    // Use distPhi to decide fate
    //
    if (distPhi < -tolerance) return VUSolid::eInside;
    if (distPhi < tolerance) return VUSolid::eSurface;
    return VUSolid::eOutside;
  }
  else
  {
    //
    // We're outside the extent of the face,
    // so the distance is penalized by distance from edges in RZ
    //
    *bestDistance = std::sqrt(distPhi * distPhi + distRZ2);

    //
    // Use edge normal to decide fate
    //
    UVector3 cc(base3Dnorm->r * radial.x(),
                base3Dnorm->r * radial.y(),
                base3Dnorm->z);
    cc = p - cc;
    double normDist = head3Dnorm->Dot(cc);
    if (distRZ2 > tolerance * tolerance)
    {
      //
      // We're far enough away that eSurface is not possible
      //
      return normDist < 0 ? VUSolid::eInside : VUSolid::eOutside;
    }

    if (normDist < -tolerance) return VUSolid::eInside;
    if (normDist <  tolerance) return VUSolid::eSurface;
    return VUSolid::eOutside;
  }
}


//
// Normal
//
// This virtual member is simple for our planer shape,
// which has only one normal
//
UVector3 UPolyPhiFace::Normal(const UVector3& p,
                              double* bestDistance)
{
  //
  // Get distance along phi, which if negative means the point
  // is nominally inside the shape.
  //
  double distPhi = normal.Dot(p);

  //
  // Calculate projected point in r,z
  //
  double r = radial.Dot(p);

  //
  // Are we inside the face?
  //
  double distRZ2;

  if (InsideEdges(r, p.z(), &distRZ2, 0))
  {
    //
    // Yup, answer is just distPhi
    //
    *bestDistance = std::fabs(distPhi);
  }
  else
  {
    //
    // Nope. Penalize by distance out
    //
    *bestDistance = std::sqrt(distPhi * distPhi + distRZ2);
  }

  return normal;
}


//
// Extent
//
// This actually isn't needed by polycone or polyhedra...
//
double UPolyPhiFace::Extent(const UVector3 axis)
{
  double max = -UUtils::kInfinity;

  UPolyPhiFaceVertex* corner = corners;
  do
  {
    double here = axis.x() * corner->r * radial.x()
                  + axis.y() * corner->r * radial.y()
                  + axis.z() * corner->z;
    if (here > max) max = here;
  }
  while (++corner < corners + numEdges);

  return max;
}


/*
//
// CalculateExtent
//
// See notes in UVCSGface
//
void UPolyPhiFace::CalculateExtent( const EAxisType axis,
                                     const UVoxelLimits &voxelLimit,
                                     const UAffineTransform &transform,
                                           USolidExtentList &extentList )
{
  //
  // Construct a (sometimes big) clippable polygon,
  //
  // Perform the necessary transformations while doing so
  //
  UClippablePolygon polygon;

  UPolyPhiFaceVertex *corner = corners;
  do
  {
    UVector3 point( 0, 0, corner->z );
    point += radial*corner->r;

    polygon.AddVertexInOrder( transform.TransformPoint( point ) );
  } while( ++corner < corners + numEdges );

  //
  // Clip away
  //
  if (polygon.PartialClip( voxelLimit, axis ))
  {
    //
    // Add it to the list
    //
    polygon.SetNormal( transform.TransformAxis(normal) );
    extentList.AddSurface( polygon );
  }
}
*/

//
//-------------------------------------------------------


//
// InsideEdgesExact
//
// Decide if the point in r,z is inside the edges of our face,
// **but** do so consistently with other faces.
//
// This routine has functionality similar to InsideEdges, but uses
// an algorithm to decide if a trajectory falls inside or outside the
// face that uses only the trajectory p,v values and the three dimensional
// points representing the edges of the polygon. The objective is to plug up
// any leaks between touching UPolyPhiFaces (at r==0) and any other face
// that uses the same convention.
//
// See: "Computational Geometry in C (Second Edition)"
// http://cs.smith.edu/~orourke/
//
bool UPolyPhiFace::InsideEdgesExact(double r, double z,
                                    double normSign,
                                    const UVector3& p,
                                    const UVector3& v)
{
  //
  // Quick check of extent
  //
  if ((r < rMin - VUSolid::Tolerance())
      || (r > rMax + VUSolid::Tolerance())) return false;

  if ((z < zMin - VUSolid::Tolerance())
      || (z > zMax + VUSolid::Tolerance())) return false;

  //
  // Exact check: loop over all vertices
  //
  double qx = p.x() + v.x(),
         qy = p.y() + v.y(),
         qz = p.z() + v.z();

  int answer = 0;
  UPolyPhiFaceVertex* corn = corners,
                      *prev = corners + numEdges - 1;

  double cornZ, prevZ;

  prevZ = ExactZOrder(z, qx, qy, qz, v, normSign, prev);
  do
  {
    //
    // Get z order of this vertex, and compare to previous vertex
    //
    cornZ = ExactZOrder(z, qx, qy, qz, v, normSign, corn);

    if (cornZ < 0)
    {
      if (prevZ < 0) continue;
    }
    else if (cornZ > 0)
    {
      if (prevZ > 0) continue;
    }
    else
    {
      //
      // By chance, we overlap exactly (within precision) with
      // the current vertex. Continue if the same happened previously
      // (e.g. the previous vertex had the same z value)
      //
      if (prevZ == 0) continue;

      //
      // Otherwise, to decide what to do, we need to know what is
      // coming up next. Specifically, we need to find the next vertex
      // with a non-zero z order.
      //
      // One might worry about infinite loops, but the above conditional
      // should prevent it
      //
      UPolyPhiFaceVertex* next = corn;
      double nextZ;
      do
      {
        next++;
        if (next == corners + numEdges) next = corners;

        nextZ = ExactZOrder(z, qx, qy, qz, v, normSign, next);
      }
      while (nextZ == 0);

      //
      // If we won't be changing direction, go to the next vertex
      //
      if (nextZ * prevZ < 0) continue;
    }


    //
    // We overlap in z with the side of the face that stretches from
    // vertex "prev" to "corn". On which side (left or right) do
    // we lay with respect to this segment?
    //
    UVector3 qa(qx - prev->x, qy - prev->y, qz - prev->z),
             qb(qx - corn->x, qy - corn->y, qz - corn->z);

    double aboveOrBelow = normSign * qa.Cross(qb).Dot(v);

    if (aboveOrBelow > 0)
      answer++;
    else if (aboveOrBelow < 0)
      answer--;
    else
    {
      //
      // A precisely zero answer here means we exactly
      // intersect (within roundoff) the edge of the face.
      // Return true in this case.
      //
      return true;
    }
  }
  while (prevZ = cornZ, prev = corn, ++corn < corners + numEdges);

//  int fanswer = std::abs(answer);
//  if (fanswer==1 || fanswer>2) {
//    cerr << "UPolyPhiFace::InsideEdgesExact: answer is "
//           << answer << std::endl;
//  }

  return answer != 0;
}


//
// InsideEdges (don't care aboud distance)
//
// Decide if the point in r,z is inside the edges of our face
//
// This routine can be made a zillion times quicker by implementing
// better code, for example:
//
//    int pnpoly(int npol, float *xp, float *yp, float x, float y)
//    {
//      int i, j, c = 0;
//      for (i = 0, j = npol-1; i < npol; j = i++) {
//        if ((((yp[i]<=y) && (y<yp[j])) ||
//             ((yp[j]<=y) && (y<yp[i]))) &&
//            (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
//
//          c = !c;
//      }
//      return c;
//    }
//
// See "Point in Polyon Strategies", Eric Haines [Graphic Gems IV]  pp. 24-46
//
// My algorithm below is rather unique, but is based on code needed to
// calculate the distance to the shape. I left it in here because ...
// well ... to test it better.
//
bool UPolyPhiFace::InsideEdges(double r, double z)
{
  //
  // Quick check of extent
  //
  if (r < rMin || r > rMax) return false;
  if (z < zMin || z > zMax) return false;

  //
  // More thorough check
  //
  double notUsed;

  return InsideEdges(r, z, &notUsed, 0);
}


//
// InsideEdges (care about distance)
//
// Decide if the point in r,z is inside the edges of our face
//
bool UPolyPhiFace::InsideEdges(double r, double z,
                               double* bestDist2,
                               UPolyPhiFaceVertex** base3Dnorm,
                               UVector3** head3Dnorm)
{
  double bestDistance2 = UUtils::kInfinity;
  bool   answer = 0;

  UPolyPhiFaceEdge* edge = edges;
  do
  {
    UPolyPhiFaceVertex* testMe;
    //
    // Get distance perpendicular to the edge
    //
    double dr = (r - edge->v0->r), dz = (z - edge->v0->z);

    double distOut = dr * edge->tz - dz * edge->tr;
    double distance2 = distOut * distOut;
    if (distance2 > bestDistance2) continue;        // No hope!

    //
    // Check to see if normal intersects edge within the edge's boundary
    //
    double q = dr * edge->tr + dz * edge->tz;

    //
    // If it doesn't, penalize distance2 appropriately
    //
    if (q < 0)
    {
      distance2 += q * q;
      testMe = edge->v0;
    }
    else if (q > edge->length)
    {
      double s2 = q - edge->length;
      distance2 += s2 * s2;
      testMe = edge->v1;
    }
    else
    {
      testMe = 0;
    }

    //
    // Closest edge so far?
    //
    if (distance2 < bestDistance2)
    {
      bestDistance2 = distance2;
      if (testMe)
      {
        double distNorm = dr * testMe->rNorm + dz * testMe->zNorm;
        answer = (distNorm <= 0);
        if (base3Dnorm)
        {
          *base3Dnorm = testMe;
          *head3Dnorm = &testMe->norm3D;
        }
      }
      else
      {
        answer = (distOut <= 0);
        if (base3Dnorm)
        {
          *base3Dnorm = edge->v0;
          *head3Dnorm = &edge->norm3D;
        }
      }
    }
  }
  while (++edge < edges + numEdges);

  *bestDist2 = bestDistance2;
  return answer;
}

//
// Calculation of Surface Area of a Triangle
// In the same time Random Point in Triangle is given
//
double UPolyPhiFace::SurfaceTriangle(UVector3 p1,
                                     UVector3 p2,
                                     UVector3 p3,
                                     UVector3* p4)
{
  UVector3 v, w;

  v = p3 - p1;
  w = p1 - p2;
  double lambda1 = UUtils::Random();
  double lambda2 = lambda1 * UUtils::Random();

  *p4 = p2 + lambda1 * w + lambda2 * v;
  return 0.5 * (v.Cross(w)).Mag();
}

//
// Compute surface area
//
double UPolyPhiFace::SurfaceArea()
{
  if (fSurfaceArea == 0.)
  {
    Triangulate();
  }
  return fSurfaceArea;
}

//
// Return random point on face
//
UVector3 UPolyPhiFace::GetPointOnFace()
{
  Triangulate();
  return surface_point;
}

//
// Auxiliary Functions used for Finding the PointOnFace using Triangulation
//

//
// Calculation of 2*Area of Triangle with Sign
//
double UPolyPhiFace::Area2(UVector2 a,
                           UVector2 b,
                           UVector2 c)
{
  return ((b.x - a.x) * (c.y - a.y) -
          (c.x - a.x) * (b.y - a.y));
}

//
// Boolean function for sign of Surface
//
bool UPolyPhiFace::Left(UVector2 a,
                        UVector2 b,
                        UVector2 c)
{
  return Area2(a, b, c) > 0;
}

//
// Boolean function for sign of Surface
//
bool UPolyPhiFace::LeftOn(UVector2 a,
                          UVector2 b,
                          UVector2 c)
{
  return Area2(a, b, c) >= 0;
}

//
// Boolean function for sign of Surface
//
bool UPolyPhiFace::Collinear(UVector2 a,
                             UVector2 b,
                             UVector2 c)
{
  return Area2(a, b, c) == 0;
}

//
// Boolean function for finding "Proper" Intersection
// That means Intersection of two lines segments (a,b) and (c,d)
//
bool UPolyPhiFace::IntersectProp(UVector2 a,
                                 UVector2 b,
                                 UVector2 c, UVector2 d)
{
  if (Collinear(a, b, c) || Collinear(a, b, d) ||
      Collinear(c, d, a) || Collinear(c, d, b))
  {
    return false;
  }

  bool Positive;
  Positive = !(Left(a, b, c)) ^ !(Left(a, b, d));
  return Positive && (!Left(c, d, a) ^ !Left(c, d, b));
}

//
// Boolean function for determining if Point c is between a and b
// For the tree points(a,b,c) on the same line
//
bool UPolyPhiFace::Between(UVector2 a, UVector2 b, UVector2 c)
{
  if (!Collinear(a, b, c))
  {
    return false;
  }

  if (a.x != b.x)
  {
    return ((a.x <= c.x) && (c.x <= b.x)) ||
           ((a.x >= c.x) && (c.x >= b.x));
  }
  else
  {
    return ((a.y <= c.y) && (c.y <= b.y)) ||
           ((a.y >= c.y) && (c.y >= b.y));
  }
}

//
// Boolean function for finding Intersection "Proper" or not
// Between two line segments (a,b) and (c,d)
//
bool UPolyPhiFace::Intersect(UVector2 a,
                             UVector2 b,
                             UVector2 c, UVector2 d)
{
  if (IntersectProp(a, b, c, d))
  {
    return true;
  }
  else if (Between(a, b, c) ||
           Between(a, b, d) ||
           Between(c, d, a) ||
           Between(c, d, b))
  {
    return true;
  }
  else
  {
    return false;
  }
}

//
// Boolean Diagonalie help to determine
// if diagonal s of segment (a,b) is convex or reflex
//
bool UPolyPhiFace::Diagonalie(UPolyPhiFaceVertex* a,
                              UPolyPhiFaceVertex* b)
{
  UPolyPhiFaceVertex*   corner = triangles;
  UPolyPhiFaceVertex*   corner_next = triangles;

  // For each Edge (corner,corner_next)
  do
  {
    corner_next = corner->next;

    // Skip edges incident to a of b
    //
    if ((corner != a) && (corner_next != a)
        && (corner != b) && (corner_next != b))
    {
      UVector2 rz1, rz2, rz3, rz4;
      rz1 = UVector2(a->r, a->z);
      rz2 = UVector2(b->r, b->z);
      rz3 = UVector2(corner->r, corner->z);
      rz4 = UVector2(corner_next->r, corner_next->z);
      if (Intersect(rz1, rz2, rz3, rz4))
      {
        return false;
      }
    }
    corner = corner->next;

  }
  while (corner != triangles);

  return true;
}

//
// Boolean function that determine if b is Inside Cone (a0,a,a1)
// being a the center of the Cone
//
bool UPolyPhiFace::InCone(UPolyPhiFaceVertex* a, UPolyPhiFaceVertex* b)
{
  // a0,a and a1 are consecutive vertices
  //
  UPolyPhiFaceVertex* a0, *a1;
  a1 = a->next;
  a0 = a->prev;

  UVector2 arz, arz0, arz1, brz;
  arz = UVector2(a->r, a->z);
  arz0 = UVector2(a0->r, a0->z);
  arz1 = UVector2(a1->r, a1->z);
  brz = UVector2(b->r, b->z);


  if (LeftOn(arz, arz1, arz0)) // If a is convex vertex
  {
    return Left(arz, brz, arz0) && Left(brz, arz, arz1);
  }
  else                       // Else a is reflex
  {
    return !(LeftOn(arz, brz, arz1) && LeftOn(brz, arz, arz0));
  }
}

//
// Boolean function finding if Diagonal is possible
// inside Polycone or PolyHedra
//
bool UPolyPhiFace::Diagonal(UPolyPhiFaceVertex* a, UPolyPhiFaceVertex* b)
{
  return InCone(a, b) && InCone(b, a) && Diagonalie(a, b);
}

//
// Initialisation for Triangulisation by ear tips
// For details see "Computational Geometry in C" by Joseph O'Rourke
//
void UPolyPhiFace::EarInit()
{
  UPolyPhiFaceVertex*   corner = triangles;
  UPolyPhiFaceVertex* c_prev, *c_next;

  do
  {
    // We need to determine three consecutive vertices
    //
    c_next = corner->next;
    c_prev = corner->prev;

    // Calculation of ears
    //
    corner->ear = Diagonal(c_prev, c_next);
    corner = corner->next;

  }
  while (corner != triangles);
}

//
// Triangulisation by ear tips for Polycone or Polyhedra
// For details see "Computational Geometry in C" by Joseph O'Rourke
//
void UPolyPhiFace::Triangulate()
{
  // The copy of Polycone is made and this copy is reordered in order to
  // have a list of triangles. This list is used for GetPointOnFace().

  UPolyPhiFaceVertex* tri_help = new UPolyPhiFaceVertex[numEdges];
  triangles = tri_help;
  UPolyPhiFaceVertex* triang = triangles;

  std::vector<double> areas;
  std::vector<UVector3> points;
  double area = 0.;
  UPolyPhiFaceVertex* v0, *v1, *v2, *v3, *v4;
  v2 = triangles;

  // Make copy for prev/next for triang=corners
  //
  UPolyPhiFaceVertex* helper = corners;
  UPolyPhiFaceVertex* helper2 = corners;
  do
  {
    triang->r = helper->r;
    triang->z = helper->z;
    triang->x = helper->x;
    triang->y = helper->y;

    // add pointer on prev corner
    //
    if (helper == corners)
    {
      triang->prev = triangles + numEdges - 1;
    }
    else
    {
      triang->prev = helper2;
    }

    // add pointer on next corner
    //
    if (helper < corners + numEdges - 1)
    {
      triang->next = triang + 1;
    }
    else
    {
      triang->next = triangles;
    }
    helper2 = triang;
    helper = helper->next;
    triang = triang->next;

  }
  while (helper != corners);

  EarInit();

  int n = numEdges;
  int i = 0;
  UVector3 p1, p2, p3, p4;
  const int max_n_loops = numEdges * 10000; // protection against infinite loop

  // Each step of outer loop removes one ear
  //
  while (n > 3) // Inner loop searches for one ear
  {
    v2 = triangles;
    do
    {
      if (v2->ear) // Ear found. Fill variables
      {
        // (v1,v3) is diagonal
        //
        v3 = v2->next;
        v4 = v3->next;
        v1 = v2->prev;
        v0 = v1->prev;

        // Calculate areas and points

        p1 = UVector3((v2)->x, (v2)->y, (v2)->z);
        p2 = UVector3((v1)->x, (v1)->y, (v1)->z);
        p3 = UVector3((v3)->x, (v3)->y, (v3)->z);

        double result1 = SurfaceTriangle(p1, p2, p3, &p4);
        points.push_back(p4);
        areas.push_back(result1);
        area = area + result1;

        // Update earity of diagonal endpoints
        //
        v1->ear = Diagonal(v0, v3);
        v3->ear = Diagonal(v1, v4);

        // Cut off the ear v2
        // Has to be done for a copy and not for real PolyPhiFace
        //
        v1->next = v3;
        v3->prev = v1;
        triangles = v3; // In case the head was v2
        n--;

        break; // out of inner loop
      }       // end if ear found

      v2 = v2->next;

    }
    while (v2 != triangles);

    i++;
    if (i >= max_n_loops)
    {
      UUtils::Exception("UPolyPhiFace::Triangulation()",
                        "GeomSolids0003", UFatalError, 1,
                        "Maximum number of steps is reached for triangulation!");
    }
  }  // end outer while loop

  if (v2->next)
  {
    // add last triangle
    //
    v2 = v2->next;
    p1 = UVector3((v2)->x, (v2)->y, (v2)->z);
    p2 = UVector3((v2->next)->x, (v2->next)->y, (v2->next)->z);
    p3 = UVector3((v2->prev)->x, (v2->prev)->y, (v2->prev)->z);
    double result1 = SurfaceTriangle(p1, p2, p3, &p4);
    points.push_back(p4);
    areas.push_back(result1);
    area = area + result1;
  }

  // Surface Area is stored
  //
  fSurfaceArea = area;

  // Second Step: choose randomly one surface
  //
  double chose = area * UUtils::Random();

  // Third Step: Get a point on choosen surface
  //
  double Achose1, Achose2;
  Achose1 = 0;
  Achose2 = 0.;
  i = 0;
  do
  {
    Achose2 += areas[i];
    if (chose >= Achose1 && chose < Achose2)
    {
      UVector3 point;
      point = points[i] ;
      surface_point = point;
      break;
    }
    i++;
    Achose1 = Achose2;
  }
  while (i < numEdges - 2);

  delete [] tri_help;
  tri_help = 0;
}
