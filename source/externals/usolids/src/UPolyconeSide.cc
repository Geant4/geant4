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
// UPolyconeSide
//
// 19.04.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include <string>
#include <cmath>
#include <sstream>
#include "UPolyconeSide.hh"
#include "UIntersectingCone.hh"
#include "VUSolid.hh"

//
// Constructor
//
// Values for r1,z1 and r2,z2 should be specified in clockwise
// order in (r,z).
//
UPolyconeSide::UPolyconeSide(const UPolyconeSideRZ* prevRZ,
                             const UPolyconeSideRZ* tail,
                             const UPolyconeSideRZ* head,
                             const UPolyconeSideRZ* nextRZ,
                             double thePhiStart,
                             double theDeltaPhi,
                             bool thePhiIsOpen,
                             bool isAllBehind)
  : ncorners(0), corners(0)
{

  tolerance = VUSolid::Tolerance();

  fSurfaceArea = 0.0;

  //
  // Record values
  //
  r[0] = tail->r;
  z[0] = tail->z;
  r[1] = head->r;
  z[1] = head->z;

  phiIsOpen = thePhiIsOpen;
  if (phiIsOpen)
  {
    deltaPhi = theDeltaPhi;
    startPhi = thePhiStart;

    //
    // Set phi values to our conventions
    //
    while (deltaPhi < 0.0) deltaPhi += 2 * UUtils::kPi;
    while (startPhi < 0.0) startPhi += 2 * UUtils::kPi;

    //
    // Calculate corner coordinates
    //
    ncorners = 4;
    corners = new UVector3[ncorners];

    corners[0] = UVector3(tail->r * std::cos(startPhi),
                          tail->r * std::sin(startPhi), tail->z);
    corners[1] = UVector3(head->r * std::cos(startPhi),
                          head->r * std::sin(startPhi), head->z);
    corners[2] = UVector3(tail->r * std::cos(startPhi + deltaPhi),
                          tail->r * std::sin(startPhi + deltaPhi), tail->z);
    corners[3] = UVector3(head->r * std::cos(startPhi + deltaPhi),
                          head->r * std::sin(startPhi + deltaPhi), head->z);
  }
  else
  {
    deltaPhi = 2 * UUtils::kPi;
    startPhi = 0.0;
  }

  allBehind = isAllBehind;

  //
  // Make our intersecting cone
  //
  cone = new UIntersectingCone(r, z);

  //
  // Calculate vectors in r,z space
  //
  rS = r[1] - r[0];
  zS = z[1] - z[0];
  length = std::sqrt(rS * rS + zS * zS);
  rS /= length;
  zS /= length;

  rNorm = +zS;
  zNorm = -rS;

  double lAdj;

  prevRS = r[0] - prevRZ->r;
  prevZS = z[0] - prevRZ->z;
  lAdj = std::sqrt(prevRS * prevRS + prevZS * prevZS);
  prevRS /= lAdj;
  prevZS /= lAdj;

  rNormEdge[0] = rNorm + prevZS;
  zNormEdge[0] = zNorm - prevRS;
  lAdj = std::sqrt(rNormEdge[0] * rNormEdge[0] + zNormEdge[0] * zNormEdge[0]);
  rNormEdge[0] /= lAdj;
  zNormEdge[0] /= lAdj;

  nextRS = nextRZ->r - r[1];
  nextZS = nextRZ->z - z[1];
  lAdj = std::sqrt(nextRS * nextRS + nextZS * nextZS);
  nextRS /= lAdj;
  nextZS /= lAdj;

  rNormEdge[1] = rNorm + nextZS;
  zNormEdge[1] = zNorm - nextRS;
  lAdj = std::sqrt(rNormEdge[1] * rNormEdge[1] + zNormEdge[1] * zNormEdge[1]);
  rNormEdge[1] /= lAdj;
  zNormEdge[1] /= lAdj;
}

//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
UPolyconeSide::UPolyconeSide(__void__&)
  : startPhi(0.), deltaPhi(0.), phiIsOpen(false), allBehind(false),
    cone(0), rNorm(0.), zNorm(0.), rS(0.), zS(0.), length(0.),
    prevRS(0.), prevZS(0.), nextRS(0.), nextZS(0.), ncorners(0), corners(0),
    tolerance(0.), fSurfaceArea(0.)
{
  r[0] = r[1] = 0.;
  z[0] = z[1] = 0.;
  rNormEdge[0] = rNormEdge[1] = 0.;
  zNormEdge[0] = zNormEdge[1] = 0.;
}

//
// Destructor
//
UPolyconeSide::~UPolyconeSide()
{
  delete cone;
  if (phiIsOpen)
  {
    delete [] corners;
  }
}


//
// Copy constructor
//
UPolyconeSide::UPolyconeSide(const UPolyconeSide& source)
  : UVCSGface(), ncorners(0), corners(0)
{

  CopyStuff(source);
}


//
// Assignment operator
//
UPolyconeSide& UPolyconeSide::operator=(const UPolyconeSide& source)
{
  if (this == &source)
  {
    return *this;
  }

  delete cone;
  if (phiIsOpen)
  {
    delete [] corners;
  }

  CopyStuff(source);

  return *this;
}


//
// CopyStuff
//
void UPolyconeSide::CopyStuff(const UPolyconeSide& source)
{
  r[0]    = source.r[0];
  r[1]    = source.r[1];
  z[0]    = source.z[0];
  z[1]    = source.z[1];

  startPhi  = source.startPhi;
  deltaPhi  = source.deltaPhi;
  phiIsOpen = source.phiIsOpen;
  allBehind = source.allBehind;

  tolerance = source.tolerance;
  fSurfaceArea = source.fSurfaceArea;

  cone    = new UIntersectingCone(*source.cone);

  rNorm   = source.rNorm;
  zNorm   = source.zNorm;
  rS    = source.rS;
  zS    = source.zS;
  length    = source.length;
  prevRS    = source.prevRS;
  prevZS    = source.prevZS;
  nextRS    = source.nextRS;
  nextZS    = source.nextZS;

  rNormEdge[0]   = source.rNormEdge[0];
  rNormEdge[1]  = source.rNormEdge[1];
  zNormEdge[0]  = source.zNormEdge[0];
  zNormEdge[1]  = source.zNormEdge[1];

  if (phiIsOpen)
  {
    ncorners = 4;
    corners = new UVector3[ncorners];

    corners[0] = source.corners[0];
    corners[1] = source.corners[1];
    corners[2] = source.corners[2];
    corners[3] = source.corners[3];
  }
}


//
// Intersect
//
bool UPolyconeSide::Distance(const UVector3& p,
                             const UVector3& v,
                             bool outgoing,
                             double surfTolerance,
                             double& distance,
                             double& distFromSurface,
                             UVector3& normal,
                             bool& isAllBehind)
{
  double s1, s2;
  double normSign = outgoing ? +1 : -1;

  isAllBehind = allBehind;

  //
  // Check for two possible intersections
  //
  int nside = cone->LineHitsCone(p, v, s1, s2);
  if (nside == 0) return false;

  //
  // Check the first side first, since it is (supposed to be) closest
  //
  UVector3 hit = p + s1 * v;

  if (PointOnCone(hit, normSign, p, v, normal))
  {
    //
    // Good intersection! What about the normal?
    //
    if (normSign * v.Dot(normal) > 0)
    {
      //
      // We have a valid intersection, but it could very easily
      // be behind the point. To decide if we tolerate this,
      // we have to see if the point p is on the surface near
      // the intersecting point.
      //
      // What does it mean exactly for the point p to be "near"
      // the intersection? It means that if we draw a line from
      // p to the hit, the line remains entirely within the
      // tolerance bounds of the cone. To test this, we can
      // ask if the normal is correct near p.
      //
      double pr = p.Perp();
      if (pr < DBL_MIN) pr = DBL_MIN;
      UVector3 pNormal(rNorm * p.x() / pr, rNorm * p.y() / pr, zNorm);
      if (normSign * v.Dot(pNormal) > 0)
      {
        //
        // p and intersection in same hemisphere
        //
        double distOutside2;
        distFromSurface = -normSign * DistanceAway(p, false, distOutside2);
        if (distOutside2 < surfTolerance * surfTolerance)
        {
          if (distFromSurface > -surfTolerance)
          {
            //
            // We are just inside or away from the
            // surface. Accept *any* value of distance.
            //
            distance = s1;
            return true;
          }
        }
      }
      else
        distFromSurface = s1;

      //
      // Accept positive distances
      //
      if (s1 > 0)
      {
        distance = s1;
        return true;
      }
    }
  }

  if (nside == 1) return false;

  //
  // Well, try the second hit
  //
  hit = p + s2 * v;

  if (PointOnCone(hit, normSign, p, v, normal))
  {
    //
    // Good intersection! What about the normal?
    //
    if (normSign * v.Dot(normal) > 0)
    {
      double pr = p.Perp();
      if (pr < DBL_MIN) pr = DBL_MIN;
      UVector3 pNormal(rNorm * p.x() / pr, rNorm * p.y() / pr, zNorm);
      if (normSign * v.Dot(pNormal) > 0)
      {
        double distOutside2;
        distFromSurface = -normSign * DistanceAway(p, false, distOutside2);
        if (distOutside2 < surfTolerance * surfTolerance)
        {
          if (distFromSurface > -surfTolerance)
          {
            distance = s2;
            return true;
          }
        }
      }
      else
        distFromSurface = s2;

      if (s2 > 0)
      {
        distance = s2;
        return true;
      }
    }
  }

  //
  // Better luck next time
  //
  return false;
}


double UPolyconeSide::Safety(const UVector3& p, bool outgoing)
{
  double normSign = outgoing ? -1 : +1;
  double distFrom, distOut2;

  //
  // We have two tries for each hemisphere. Try the closest first.
  //
  distFrom = normSign * DistanceAway(p, false, distOut2);
  if (distFrom > -0.5 * VUSolid::Tolerance())
  {
    //
    // Good answer
    //
    if (distOut2 > 0)
      return std::sqrt(distFrom * distFrom + distOut2);
    else
      return std::fabs(distFrom);
  }

  //
  // Try second side.
  //
  distFrom = normSign * DistanceAway(p, true, distOut2);
  if (distFrom > -0.5 * VUSolid::Tolerance())
  {

    if (distOut2 > 0)
      return std::sqrt(distFrom * distFrom + distOut2);
    else
      return std::fabs(distFrom);
  }

  return UUtils::kInfinity;
}


//
// Inside
//
VUSolid::EnumInside UPolyconeSide::Inside(const UVector3& p,
                                          double atolerance,
                                          double* bestDistance)
{
  //
  // Check both sides
  //
  double distFrom, distOut2, dist2;
  double edgeRZnorm;

  distFrom = DistanceAway(p, distOut2, &edgeRZnorm);
  dist2 = distFrom * distFrom + distOut2;


  *bestDistance = std::sqrt(dist2);   // could sqrt be removed?

  //
  // Okay then, inside or out?
  //
  if ((std::fabs(edgeRZnorm) < atolerance)
      && (distOut2 < atolerance * atolerance))
    return VUSolid::eSurface;
  else if (edgeRZnorm < 0)
    return VUSolid::eInside;
  else
    return VUSolid::eOutside;
}


//
// Normal
//
UVector3 UPolyconeSide::Normal(const UVector3& p,
                               double* bestDistance)
{
  if (p == UVector3(0., 0., 0.))
  {
    return p;
  }

  double dFrom, dOut2;

  dFrom = DistanceAway(p, false, dOut2);

  *bestDistance = std::sqrt(dFrom * dFrom + dOut2);

  double rds = p.Perp();
  if (rds != 0.)
  {
    return UVector3(rNorm * p.x() / rds, rNorm * p.y() / rds, zNorm);
  }
  return UVector3(0., 0., zNorm).Unit();
}


//
// Extent
//
double UPolyconeSide::Extent(const UVector3 axis)
{
  if (axis.Perp2() < DBL_MIN)
  {
    //
    // Special case
    //
    return axis.z() < 0 ? -cone->ZLo() : cone->ZHi();
  }

  //
  // Is the axis pointing inside our phi gap?
  //
  if (phiIsOpen)
  {
    double phi = GetPhi(axis);
    while (phi < startPhi) phi += 2 * UUtils::kPi;

    if (phi > deltaPhi + startPhi)
    {
      //
      // Yeah, looks so. Make four three vectors defining the phi
      // opening
      //
      double cosP = std::cos(startPhi), sinP = std::sin(startPhi);
      UVector3 a(r[0]*cosP, r[0]*sinP, z[0]);
      UVector3 b(r[1]*cosP, r[1]*sinP, z[1]);
      cosP = std::cos(startPhi + deltaPhi);
      sinP = std::sin(startPhi + deltaPhi);
      UVector3 c(r[0]*cosP, r[0]*sinP, z[0]);
      UVector3 d(r[1]*cosP, r[1]*sinP, z[1]);

      double ad = axis.Dot(a),
             bd = axis.Dot(b),
             cd = axis.Dot(c),
             dd = axis.Dot(d);

      if (bd > ad) ad = bd;
      if (cd > ad) ad = cd;
      if (dd > ad) ad = dd;

      return ad;
    }
  }

  //
  // Check either end
  //
  double aPerp = axis.Perp();

  double a = aPerp * r[0] + axis.z() * z[0];
  double b = aPerp * r[1] + axis.z() * z[1];

  if (b > a) a = b;

  return a;
}



//
// CalculateExtent
//
// See notes in UVCSGface
//

/*
void UPolyconeSide::CalculateExtent( const EAxisType axis,
                                      const UVoxelLimits &voxelLimit,
                                      const UAffineTransform &transform,
                                            USolidExtentList &extentList )
{
  UClippablePolygon polygon;

  //
  // Here we will approximate (ala UCons) and divide our conical section
  // into segments, like UPolyhedra. When doing so, the radius
  // is extented far enough such that the segments always lie
  // just outside the surface of the conical section we are
  // approximating.
  //

  //
  // Choose phi size of our segment(s) based on constants as
  // defined in meshdefs.hh
  //
  int numPhi = (int)(deltaPhi/UUtils::kMeshAngleDefault) + 1;
  if (numPhi < UUtils::kMinMeshSections)
    numPhi = UUtils::kMinMeshSections;
  else if (numPhi > UUtils::kMaxMeshSections)
    numPhi = UUtils::kMaxMeshSections;

  double sigPhi = deltaPhi/numPhi;

  //
  // Determine radius factor to keep segments outside
  //
  double rFudge = 1.0/std::cos(0.5*sigPhi);

  //
  // Decide which radius to use on each end of the side,
  // and whether a transition mesh is required
  //
  // {r0,z0}  - Beginning of this side
  // {r1,z1}  - Ending of this side
  // {r2,z0}  - Beginning of transition piece connecting previous
  //            side (and ends at beginning of this side)
  //
  // So, order is 2 --> 0 --> 1.
  //                    -------
  //
  // r2 < 0 indicates that no transition piece is required
  //
  double r0, r1, r2, z0, z1;

  r2 = -1;  // By default: no transition piece

  if (rNorm < -DBL_MIN)
  {
    //
    // This side faces *inward*, and so our mesh has
    // the same radius
    //
    r1 = r[1];
    z1 = z[1];
    z0 = z[0];
    r0 = r[0];

    r2 = -1;

    if (prevZS > DBL_MIN)
    {
      //
      // The previous side is facing outwards
      //
      if ( prevRS*zS - prevZS*rS > 0 )
      {
        //
        // Transition was convex: build transition piece
        //
        if (r[0] > DBL_MIN) r2 = r[0]*rFudge;
      }
      else
      {
        //
        // Transition was concave: short this side
        //
        FindLineIntersect( z0, r0, zS, rS,
                           z0, r0*rFudge, prevZS, prevRS*rFudge, z0, r0 );
      }
    }

    if ( nextZS > DBL_MIN && (rS*nextZS - zS*nextRS < 0) )
    {
      //
      // The next side is facing outwards, forming a
      // concave transition: short this side
      //
      FindLineIntersect( z1, r1, zS, rS,
                         z1, r1*rFudge, nextZS, nextRS*rFudge, z1, r1 );
    }
  }
  else if (rNorm > DBL_MIN)
  {
    //
    // This side faces *outward* and is given a boost to
    // it radius
    //
    r0 = r[0]*rFudge;
    z0 = z[0];
    r1 = r[1]*rFudge;
    z1 = z[1];

    if (prevZS < -DBL_MIN)
    {
      //
      // The previous side is facing inwards
      //
      if ( prevRS*zS - prevZS*rS > 0 )
      {
        //
        // Transition was convex: build transition piece
        //
        if (r[0] > DBL_MIN) r2 = r[0];
      }
      else
      {
        //
        // Transition was concave: short this side
        //
        FindLineIntersect( z0, r0, zS, rS*rFudge,
                           z0, r[0], prevZS, prevRS, z0, r0 );
      }
    }

    if ( nextZS < -DBL_MIN && (rS*nextZS - zS*nextRS < 0) )
    {
      //
      // The next side is facing inwards, forming a
      // concave transition: short this side
      //
      FindLineIntersect( z1, r1, zS, rS*rFudge,
                         z1, r[1], nextZS, nextRS, z1, r1 );
    }
  }
  else
  {
    //
    // This side is perpendicular to the z axis (is a disk)
    //
    // Whether or not r0 needs a rFudge factor depends
    // on the normal of the previous edge. Similar with r1
    // and the next edge. No transition piece is required.
    //
    r0 = r[0];
    r1 = r[1];
    z0 = z[0];
    z1 = z[1];

    if (prevZS > DBL_MIN) r0 *= rFudge;
    if (nextZS > DBL_MIN) r1 *= rFudge;
  }

  //
  // Loop
  //
  double phi = startPhi,
           cosPhi = std::cos(phi),
           sinPhi = std::sin(phi);

  UVector3 v0( r0*cosPhi, r0*sinPhi, z0 ),
                    v1( r1*cosPhi, r1*sinPhi, z1 ),
  v2, w0, w1, w2;
  transform.ApplyPointTransform( v0 );
  transform.ApplyPointTransform( v1 );

  if (r2 >= 0)
  {
    v2 = UVector3( r2*cosPhi, r2*sinPhi, z0 );
    transform.ApplyPointTransform( v2 );
  }

  do
  {
    phi += sigPhi;
    if (numPhi == 1) phi = startPhi+deltaPhi; // Try to avoid roundoff
    cosPhi = std::cos(phi),
    sinPhi = std::sin(phi);

    w0 = UVector3( r0*cosPhi, r0*sinPhi, z0 );
    w1 = UVector3( r1*cosPhi, r1*sinPhi, z1 );
    transform.ApplyPointTransform( w0 );
    transform.ApplyPointTransform( w1 );

    UVector3 deltaV = r0 > r1 ? w0-v0 : w1-v1;

    //
    // Build polygon, taking special care to keep the vertices
    // in order
    //
    polygon.ClearAllVertices();

    polygon.AddVertexInOrder( v0 );
    polygon.AddVertexInOrder( v1 );
    polygon.AddVertexInOrder( w1 );
    polygon.AddVertexInOrder( w0 );

    //
    // Get extent
    //
    if (polygon.PartialClip( voxelLimit, axis ))
    {
      //
      // Get Dot product of normal with target axis
      //
      polygon.SetNormal( deltaV.Cross(v1-v0).Unit() );

      extentList.AddSurface( polygon );
    }

    if (r2 >= 0)
    {
      //
      // Repeat, for transition piece
      //
      w2 = UVector3( r2*cosPhi, r2*sinPhi, z0 );
      transform.ApplyPointTransform( w2 );

      polygon.ClearAllVertices();

      polygon.AddVertexInOrder( v2 );
      polygon.AddVertexInOrder( v0 );
      polygon.AddVertexInOrder( w0 );
      polygon.AddVertexInOrder( w2 );

      if (polygon.PartialClip( voxelLimit, axis ))
      {
        polygon.SetNormal( deltaV.Cross(v0-v2).Unit() );

        extentList.AddSurface( polygon );
      }

      v2 = w2;
    }

    //
    // Next vertex
    //
    v0 = w0;
    v1 = w1;
  } while( --numPhi > 0 );

  //
  // We are almost done. But, it is important that we leave no
  // gaps in the surface of our solid. By using rFudge, however,
  // we've done exactly that, if we have a phi segment.
  // Add two additional faces if necessary
  //
  if (phiIsOpen && rNorm > DBL_MIN)
  {
    cosPhi = std::cos(startPhi);
    sinPhi = std::sin(startPhi);

    UVector3 a0( r[0]*cosPhi, r[0]*sinPhi, z[0] ),
                  a1( r[1]*cosPhi, r[1]*sinPhi, z[1] ),
                  b0( r0*cosPhi, r0*sinPhi, z[0] ),
                  b1( r1*cosPhi, r1*sinPhi, z[1] );

    transform.ApplyPointTransform( a0 );
    transform.ApplyPointTransform( a1 );
    transform.ApplyPointTransform( b0 );
    transform.ApplyPointTransform( b1 );

    polygon.ClearAllVertices();

    polygon.AddVertexInOrder( a0 );
    polygon.AddVertexInOrder( a1 );
    polygon.AddVertexInOrder( b0 );
    polygon.AddVertexInOrder( b1 );

    if (polygon.PartialClip( voxelLimit , axis))
    {
      UVector3 normal( sinPhi, -cosPhi, 0 );
      polygon.SetNormal( transform.TransformAxis( normal ) );

      extentList.AddSurface( polygon );
    }

    cosPhi = std::cos(startPhi+deltaPhi);
    sinPhi = std::sin(startPhi+deltaPhi);

    a0 = UVector3( r[0]*cosPhi, r[0]*sinPhi, z[0] ),
    a1 = UVector3( r[1]*cosPhi, r[1]*sinPhi, z[1] ),
    b0 = UVector3( r0*cosPhi, r0*sinPhi, z[0] ),
    b1 = UVector3( r1*cosPhi, r1*sinPhi, z[1] );
    transform.ApplyPointTransform( a0 );
    transform.ApplyPointTransform( a1 );
    transform.ApplyPointTransform( b0 );
    transform.ApplyPointTransform( b1 );

    polygon.ClearAllVertices();

    polygon.AddVertexInOrder( a0 );
    polygon.AddVertexInOrder( a1 );
    polygon.AddVertexInOrder( b0 );
    polygon.AddVertexInOrder( b1 );

    if (polygon.PartialClip( voxelLimit, axis ))
    {
      UVector3 normal( -sinPhi, cosPhi, 0 );
      polygon.SetNormal( transform.TransformAxis( normal ) );

      extentList.AddSurface( polygon );
    }
  }

  return;
}
*/

//
// GetPhi
//
// Calculate Phi for a given 3-vector (point), if not already cached for the
// same point, in the attempt to avoid consecutive computation of the same
// quantity
//
double UPolyconeSide::GetPhi(const UVector3& p)
{
  double val = 0.;

  val = p.Phi();

  return val;
}

//
// DistanceAway
//
// Calculate distance of a point from our conical surface, including the effect
// of any phi segmentation
//
// Arguments:
//  p            - (in) Point to check
//  opposite      - (in) If true, check opposite hemisphere (see below)
//  distOutside  - (out) Additional distance outside the edges of the surface
//  edgeRZnorm    - (out) if negative, point is inside
//
//  return value = distance from the conical plane, if extrapolated beyond edges,
//                 signed by whether the point is in inside or outside the shape
//
// Notes:
//  * There are two answers, depending on which hemisphere is considered.
//
double UPolyconeSide::DistanceAway(const UVector3& p,
                                   bool opposite,
                                   double& distOutside2,
                                   double* edgeRZnorm)
{
  //
  // Convert our point to r and z
  //
  double rx = p.Perp(), zx = p.z();

  //
  // Change sign of r if opposite says we should
  //
  if (opposite) rx = -rx;

  //
  // Calculate return value
  //
  double deltaR = rx - r[0], deltaZ = zx - z[0];
  double answer = deltaR * rNorm + deltaZ * zNorm;

  //
  // Are we off the surface in r,z space?
  //
  double q = deltaR * rS + deltaZ * zS;
  if (q < 0)
  {
    distOutside2 = q * q;
    if (edgeRZnorm) *edgeRZnorm = deltaR * rNormEdge[0] + deltaZ * zNormEdge[0];
  }
  else if (q > length)
  {
    distOutside2 = UUtils::sqr(q - length);
    if (edgeRZnorm)
    {
      deltaR = rx - r[1];
      deltaZ = zx - z[1];
      *edgeRZnorm = deltaR * rNormEdge[1] + deltaZ * zNormEdge[1];
    }
  }
  else
  {
    distOutside2 = 0;
    if (edgeRZnorm) *edgeRZnorm = answer;
  }

  if (phiIsOpen)
  {
    //
    // Finally, check phi
    //
    double phi = GetPhi(p);
    while (phi < startPhi) phi += 2 * UUtils::kPi;

    if (phi > startPhi + deltaPhi)
    {
      //
      // Oops. Are we closer to the start phi or end phi?
      //
      double d1 = phi - startPhi - deltaPhi;
      while (phi > startPhi) phi -= 2 * UUtils::kPi;
      double d2 = startPhi - phi;

      if (d2 < d1) d1 = d2;

      //
      // Add result to our distance
      //
      double dist = d1 * rx;

      distOutside2 += dist * dist;
      if (edgeRZnorm)
      {
        *edgeRZnorm = std::max(std::fabs(*edgeRZnorm), std::fabs(dist));
      }
    }
  }

  return answer;
}

//
// DistanceAway
//
//
// Special version of DistanceAway for Inside.
// opposite parameter is not used, instead use sign of rx for choosing the side
//
double UPolyconeSide::DistanceAway(const UVector3& p,
                                   double& distOutside2,
                                   double* edgeRZnorm)
{
  //
  // Convert our point to r and z
  //
  double rx = p.Perp(), zx = p.z();

  //
  // Change sign of r if we should
  //
  int part = 1;
  if (rx < 0) part = -1;

  //
  // Calculate return value
  //
  double deltaR = rx - r[0]*part, deltaZ = zx - z[0];
  double answer = deltaR * rNorm*part + deltaZ * zNorm;

  //
  // Are we off the surface in r,z space?
  //
  double q = deltaR * rS *part+ deltaZ * zS;
  if (q < 0)
  {
    distOutside2 = q * q;
    if (edgeRZnorm)
    {
      *edgeRZnorm = deltaR * rNormEdge[0]*part + deltaZ * zNormEdge[0];
    }
  }
  else if (q > length)
  {
    distOutside2 = UUtils::sqr(q - length);
    if (edgeRZnorm)
    {
      deltaR = rx - r[1]*part;
      deltaZ = zx - z[1];
      *edgeRZnorm = deltaR * rNormEdge[1]*part + deltaZ * zNormEdge[1];
    }
  }
  else
  {
    distOutside2 = 0;
    if (edgeRZnorm) *edgeRZnorm = answer;
  }

  if (phiIsOpen)
  {
    //
    // Finally, check phi
    //
    double phi = GetPhi(p);
    while (phi < startPhi) phi += 2 * UUtils::kPi;

    if (phi > startPhi + deltaPhi)
    {
      //
      // Oops. Are we closer to the start phi or end phi?
      //
      double d1 = phi - startPhi - deltaPhi;
      while (phi > startPhi) phi -= 2 * UUtils::kPi;
      double d2 = startPhi - phi;

      if (d2 < d1) d1 = d2;

      //
      // Add result to our distance
      //
      double dist = d1 * rx*part;

      distOutside2 += dist * dist;
      if (edgeRZnorm)
      {
        *edgeRZnorm = std::max(std::fabs(*edgeRZnorm), std::fabs(dist));
      }
    }
  }

  return answer;
}


//
// PointOnCone
//
// Decide if a point is on a cone and return normal if it is
//
bool UPolyconeSide::PointOnCone(const UVector3& hit,
                                double normSign,
                                const UVector3& p,
                                const UVector3& v,
                                UVector3& normal)
{
  double rx = hit.Perp();
  //
  // Check radial/z extent, as appropriate
  //
  if (!cone->HitOn(rx, hit.z())) return false;

  if (phiIsOpen)
  {
    double phiTolerant = 2.0 * VUSolid::Tolerance() / (rx + VUSolid::Tolerance());
    //
    // Check phi segment. Here we have to be careful
    // to use the standard method consistent with
    // PolyPhiFace. See PolyPhiFace::InsideEdgesExact
    //
    double phi = GetPhi(hit);
    while (phi < startPhi - phiTolerant) phi += 2 * UUtils::kPi;

    if (phi > startPhi + deltaPhi + phiTolerant) return false;

    if (phi > startPhi + deltaPhi - phiTolerant)
    {
      //
      // Exact treatment
      //
      UVector3 qx = p + v;
      UVector3 qa = qx - corners[2],
               qb = qx - corners[3];
      UVector3 qacb = qa.Cross(qb);

      if (normSign * qacb.Dot(v) < 0) return false;
    }
    else if (phi < phiTolerant)
    {
      UVector3 qx = p + v;
      UVector3 qa = qx - corners[1],
               qb = qx - corners[0];
      UVector3 qacb = qa.Cross(qb);

      if (normSign * qacb.Dot(v) < 0) return false;
    }
  }

  //
  // We have a good hit! Calculate normal
  //
  if (rx < DBL_MIN)
    normal = UVector3(0, 0, zNorm < 0 ? -1 : 1);
  else
    normal = UVector3(rNorm * hit.x() / rx, rNorm * hit.y() / rx, zNorm);
  return true;
}


//
// FindLineIntersect
//
// Decide the point at which two 2-dimensional lines intersect
//
// Equation of line: x = x1 + s*tx1
//                   y = y1 + s*ty1
//
// It is assumed that the lines are *not* parallel
//
void UPolyconeSide::FindLineIntersect(double x1,  double y1,
                                      double tx1, double ty1,
                                      double x2,  double y2,
                                      double tx2, double ty2,
                                      double& x,  double& y)
{
  //
  // The solution is a simple linear equation
  //
  double deter = tx1 * ty2 - tx2 * ty1;

  double s1 = ((x2 - x1) * ty2 - tx2 * (y2 - y1)) / deter;
  double s2 = ((x2 - x1) * ty1 - tx1 * (y2 - y1)) / deter;

  //
  // We want the answer to not depend on which order the
  // lines were specified. Take average.
  //
  x = 0.5 * (x1 + s1 * tx1 + x2 + s2 * tx2);
  y = 0.5 * (y1 + s1 * ty1 + y2 + s2 * ty2);
}

//
// Calculate surface area for GetPointOnSurface()
//
double UPolyconeSide::SurfaceArea()
{
  if (fSurfaceArea == 0)
  {
    fSurfaceArea = (r[0] + r[1]) * std::sqrt(UUtils::sqr(r[0] - r[1]) + UUtils::sqr(z[0] - z[1]));
    fSurfaceArea *= 0.5 * (deltaPhi);
  }
  return fSurfaceArea;
}

//
// GetPointOnFace
//
UVector3 UPolyconeSide::GetPointOnFace()
{
  double x, y, zz;
  double rr, phi, dz, dr;
  dr = r[1] - r[0];
  dz = z[1] - z[0];
  phi = startPhi + deltaPhi * UUtils::Random();
  rr = r[0] + dr * UUtils::Random();

  x = rr * std::cos(phi);
  y = rr * std::sin(phi);

  // PolyconeSide has a Ring Form
  //
  if (dz == 0.)
  {
    zz = z[0];
  }
  else
  {
    if (dr == 0.) // PolyconeSide has a Tube Form
    {
      zz = z[0] + dz * UUtils::Random();
    }
    else
    {
      zz = z[0] + (rr - r[0]) * dz / dr;
    }
  }

  return UVector3(x, y, zz);
}
