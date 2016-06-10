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
// UPolyhedraSide
//
// 19.09.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include <string>
#include <cmath>
#include <sstream>
#include "UPolyhedraSide.hh"
#include "UIntersectingCone.hh"


//
// Constructor
//
// Values for r1,z1 and r2,z2 should be specified in clockwise
// order in (r,z).
//
UPolyhedraSide::UPolyhedraSide(const UPolyhedraSideRZ* prevRZ,
                               const UPolyhedraSideRZ* tail,
                               const UPolyhedraSideRZ* head,
                               const UPolyhedraSideRZ* nextRZ,
                               int theNumSide,
                               double thePhiStart,
                               double thePhiTotal,
                               bool thePhiIsOpen,
                               bool isAllBehind)
{
  kCarTolerance = VUSolid::Tolerance();
  fSurfaceArea = 0.;
  fPhi.first.Set(0);
  fPhi.second = 0.0;

  //
  // Record values
  //
  r[0] = tail->r;
  z[0] = tail->z;
  r[1] = head->r;
  z[1] = head->z;

  double phiTotal;

  //
  // Set phi to our convention
  //
  startPhi = thePhiStart;
  while (startPhi < 0.0) startPhi += 2 * UUtils::kPi;

  phiIsOpen = thePhiIsOpen;
  phiTotal = (phiIsOpen) ? thePhiTotal : 2 * UUtils::kPi;

  allBehind = isAllBehind;

  //
  // Make our intersecting cone
  //
  cone = new UIntersectingCone(r, z);

  //
  // Construct side plane vector Set
  //
  numSide = theNumSide;
  deltaPhi = phiTotal / theNumSide;
  endPhi = startPhi + phiTotal;

  vecs = new UPolyhedraSideVec[numSide];

  edges = new UPolyhedraSideEdge[phiIsOpen ? numSide + 1 : numSide];

  //
  // ...this is where we start
  //
  double phi = startPhi;
  UVector3 a1(r[0]*std::cos(phi), r[0]*std::sin(phi), z[0]),
           b1(r[1]*std::cos(phi), r[1]*std::sin(phi), z[1]),
           c1(prevRZ->r * std::cos(phi), prevRZ->r * std::sin(phi), prevRZ->z),
           d1(nextRZ->r * std::cos(phi), nextRZ->r * std::sin(phi), nextRZ->z),
           a2, b2, c2, d2;
  UPolyhedraSideEdge* edge = edges;

  UPolyhedraSideVec* vec = vecs;
  do
  {
    //
    // ...this is where we are going
    //
    phi += deltaPhi;
    a2 = UVector3(r[0] * std::cos(phi), r[0] * std::sin(phi), z[0]);
    b2 = UVector3(r[1] * std::cos(phi), r[1] * std::sin(phi), z[1]);
    c2 = UVector3(prevRZ->r * std::cos(phi), prevRZ->r * std::sin(phi), prevRZ->z);
    d2 = UVector3(nextRZ->r * std::cos(phi), nextRZ->r * std::sin(phi), nextRZ->z);

    UVector3 tt;

    //
    // ...build some relevant vectors.
    //    the point is to sacrifice a little memory with precalcs
    //    to gain speed
    //
    vec->center = 0.25 * (a1 + a2 + b1 + b2);

    tt = b2 + b1 - a2 - a1;
    vec->surfRZ = tt.Unit();
    if (vec == vecs) lenRZ = 0.25 * tt.Mag();

    tt = b2 - b1 + a2 - a1;
    vec->surfPhi = tt.Unit();
    if (vec == vecs)
    {
      lenPhi[0] = 0.25 * tt.Mag();
      tt = b2 - b1;
      lenPhi[1] = (0.5 * tt.Mag() - lenPhi[0]) / lenRZ;
    }

    tt = vec->surfPhi.Cross(vec->surfRZ);
    vec->normal = tt.Unit();

    //
    // ...edge normals are the average of the normals of
    //    the two faces they connect.
    //
    // ...edge normals are necessary if we are to accurately
    //    decide if a point is "inside" a face. For non-convex
    //    shapes, it is absolutely necessary to know information
    //    on adjacent faces to accurate determine this.
    //
    // ...we don't need them for the phi edges, since that
    //    information is taken care of internally. The r/z edges,
    //    however, depend on the adjacent UPolyhedraSide.
    //
    UVector3 a12, adj;

    a12 = a2 - a1;

    adj = 0.5 * (c1 + c2 - a1 - a2);
    adj = adj.Cross(a12);
    adj = adj.Unit() + vec->normal;
    vec->edgeNorm[0] = adj.Unit();

    a12 = b1 - b2;
    adj = 0.5 * (d1 + d2 - b1 - b2);
    adj = adj.Cross(a12);
    adj = adj.Unit() + vec->normal;
    vec->edgeNorm[1] = adj.Unit();

    //
    // ...the corners are crucial. It is important that
    //    they are calculated consistently for adjacent
    //    UPolyhedraSides, to avoid gaps caused by roundoff.
    //
    vec->edges[0] = edge;
    edge->corner[0] = a1;
    edge->corner[1] = b1;
    edge++;
    vec->edges[1] = edge;

    a1 = a2;
    b1 = b2;
    c1 = c2;
    d1 = d2;
  }
  while (++vec < vecs + numSide);

  //
  // Clean up hanging edge
  //
  if (phiIsOpen)
  {
    edge->corner[0] = a2;
    edge->corner[1] = b2;
  }
  else
  {
    vecs[numSide - 1].edges[1] = edges;
  }

  //
  // Go back and fill in remaining fields in edges
  //
  vec = vecs;
  UPolyhedraSideVec* prev = vecs + numSide - 1;
  do
  {
    edge = vec->edges[0];    // The edge between prev and vec

    //
    // Okay: edge normal is average of normals of adjacent faces
    //
    UVector3 eNorm = vec->normal + prev->normal;
    edge->normal = eNorm.Unit();

    //
    // Vertex normal is average of norms of adjacent surfaces (all four)
    // However, vec->edgeNorm is Unit vector in some direction
    // as the sum of normals of adjacent PolyhedraSide with vec.
    // The normalization used for this vector should be the same
    // for vec and prev.
    //
    eNorm = vec->edgeNorm[0] + prev->edgeNorm[0];
    edge->cornNorm[0] = eNorm.Unit();

    eNorm = vec->edgeNorm[1] + prev->edgeNorm[1];
    edge->cornNorm[1] = eNorm.Unit();
  }
  while (prev = vec, ++vec < vecs + numSide);

  if (phiIsOpen)
  {
    // double rFact = std::cos(0.5*deltaPhi);
    //
    // If phi is open, we need to patch up normals of the
    // first and last edges and their corresponding
    // vertices.
    //
    // We use vectors that are in the plane of the
    // face. This should be safe.
    //
    vec = vecs;

    UVector3 normvec = vec->edges[0]->corner[0]
                       - vec->edges[0]->corner[1];
    normvec = normvec.Cross(vec->normal);
    if (normvec.Dot(vec->surfPhi) > 0) normvec = -normvec;

    vec->edges[0]->normal = normvec.Unit();

    vec->edges[0]->cornNorm[0] = (vec->edges[0]->corner[0]
                                  - vec->center).Unit();
    vec->edges[0]->cornNorm[1] = (vec->edges[0]->corner[1]
                                  - vec->center).Unit();

    //
    // Repeat for ending phi
    //
    vec = vecs + numSide - 1;

    normvec = vec->edges[1]->corner[0] - vec->edges[1]->corner[1];
    normvec = normvec.Cross(vec->normal);
    if (normvec.Dot(vec->surfPhi) < 0) normvec = -normvec;

    vec->edges[1]->normal = normvec.Unit();

    vec->edges[1]->cornNorm[0] = (vec->edges[1]->corner[0]
                                  - vec->center).Unit();
    vec->edges[1]->cornNorm[1] = (vec->edges[1]->corner[1]
                                  - vec->center).Unit();
  }

  //
  // edgeNorm is the factor one multiplies the distance along vector phi
  // on the surface of one of our sides in order to calculate the distance
  // from the edge. (see routine DistanceAway)
  //
  edgeNorm = 1.0 / std::sqrt(1.0 + lenPhi[1] * lenPhi[1]);
}


//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
UPolyhedraSide::UPolyhedraSide(__void__&)
  : numSide(0), startPhi(0.), deltaPhi(0.), endPhi(0.),
    phiIsOpen(false), allBehind(false), cone(0), vecs(0), edges(0),
    lenRZ(0.), edgeNorm(0.), kCarTolerance(0.), fSurfaceArea(0.)
{
  r[0] = r[1] = 0.;
  z[0] = z[1] = 0.;
  lenPhi[0] = lenPhi[1] = 0.;
}


//
// Destructor
//
UPolyhedraSide::~UPolyhedraSide()
{
  delete cone;
  delete [] vecs;
  delete [] edges;
}


//
// Copy constructor
//
UPolyhedraSide::UPolyhedraSide(const UPolyhedraSide& source)
  : UVCSGface()
{
  CopyStuff(source);
}


//
// Assignment operator
//
UPolyhedraSide& UPolyhedraSide::operator=(const UPolyhedraSide& source)
{
  if (this == &source) return *this;

  delete cone;
  delete [] vecs;
  delete [] edges;

  CopyStuff(source);

  return *this;
}


//
// CopyStuff
//
void UPolyhedraSide::CopyStuff(const UPolyhedraSide& source)
{
  //
  // The simple stuff
  //
  numSide    = source.numSide;
  r[0]    = source.r[0];
  r[1]    = source.r[1];
  z[0]    = source.z[0];
  z[1]    = source.z[1];
  startPhi  = source.startPhi;
  deltaPhi  = source.deltaPhi;
  endPhi    = source.endPhi;
  phiIsOpen = source.phiIsOpen;
  allBehind = source.allBehind;

  lenRZ     = source.lenRZ;
  lenPhi[0] = source.lenPhi[0];
  lenPhi[1] = source.lenPhi[1];
  edgeNorm  = source.edgeNorm;

  kCarTolerance = source.kCarTolerance;
  fSurfaceArea = source.fSurfaceArea;

  cone = new UIntersectingCone(*source.cone);

  //
  // Duplicate edges
  //
  int  numEdges = phiIsOpen ? numSide + 1 : numSide;
  edges = new UPolyhedraSideEdge[numEdges];

  UPolyhedraSideEdge* edge = edges,
                      *sourceEdge = source.edges;
  do
  {
    *edge = *sourceEdge;
  }
  while (++sourceEdge, ++edge < edges + numEdges);

  //
  // Duplicate vecs
  //
  vecs = new UPolyhedraSideVec[numSide];

  UPolyhedraSideVec* vec = vecs,
                     *sourceVec = source.vecs;
  do
  {
    *vec = *sourceVec;
    vec->edges[0] = edges + (sourceVec->edges[0] - source.edges);
    vec->edges[1] = edges + (sourceVec->edges[1] - source.edges);
  }
  while (++sourceVec, ++vec < vecs + numSide);
}


//
// Intersect
//
// Decide if a line intersects the face.
//
// Arguments:
//  p    = (in) starting point of line segment
//  v    = (in) direction of line segment (assumed a Unit vector)
//  A, B    = (in) 2d transform variables (see note top of file)
//  normSign  = (in) desired sign for Dot product with normal (see below)
//  surfTolerance  = (in) minimum distance from the surface
//  vecs    = (in) Vector Set array
//  distance  = (out) distance to surface furfilling all requirements
//  distFromSurface = (out) distance from the surface
//  thisNormal  = (out) normal vector of the intersecting surface
//
// Return value:
//  true if an intersection is found. Otherwise, output parameters are
//  undefined.
//
// Notes:
// * normSign: if we are "inside" the shape and only want to find out how far
//   to leave the shape, we only want to consider intersections with surfaces in
//   which the trajectory is leaving the shape. Since the normal vectors to the
//   surface always point outwards from the inside, this means we want the Dot
//   product of the trajectory direction v and the normal of the side normals[i]
//   to be positive. Thus, we should specify normSign as +1.0. Otherwise, if
//   we are outside and want to go in, normSign should be Set to -1.0.
//   Don't Set normSign to zero, or you will get no intersections!
//
// * surfTolerance: see notes on argument "surfTolerance" in routine
//   "IntersectSidePlane".
//   ----HOWEVER---- We should *not* apply this surface tolerance if the
//   starting point is not within phi or z of the surface. Specifically,
//   if the starting point p angle in x/y places it on a separate side from the
//   intersection or if the starting point p is outside the z bounds of the
//   segment, surfTolerance must be ignored or we should *always* accept the
//   intersection!
//   This is simply because the sides do not have infinite extent.
//
//
bool UPolyhedraSide::Distance(const UVector3& p,
                              const UVector3& v,
                              bool outgoing,
                              double surfTolerance,
                              double& distance,
                              double& distFromSurface,
                              UVector3& normal,
                              bool& isAllBehind)
{
  int segment = -1;

  //
  // Testing the intersection of individual phi faces is
  // pretty straight forward. The simple thing therefore is to
  // form a loop and check them all in sequence.
  //
  // We try to quickly decide
  // which face would be intersected. One can make a very
  // good guess by using the intersection with a cone.
  // However, this is only reliable in 99% of the cases.
  //
  // We make a decent guess as to the one or
  // two potential faces might get intersected, and then
  // test them. If we have the wrong face, use the test
  // to make a better guess.
  //
  // Since we might have two guesses, form a queue of
  // potential intersecting faces. Keep an array of
  // already tested faces to avoid doing one more than
  // once.
  //
  // Result: at worst, an iterative search. On average,
  // a little more than two tests would be required.
  //

  if (numSide > 5)
  {
    //todo: maybe we could even use the second solution? how much also the second would be relevant???
    double s1, s2;
    int solutions = cone->LineHitsCone(p, v, s1, s2);
    if (!solutions) return false;
    if (solutions == 2 && s2 > 0 && (s2 < s1 || s1 < 0))
      s1 = s2;

    segment = PhiSegment(std::atan2(p.y + s1 * v.y, p.x + s1 * v.x));
  }

  UVector3 q = p + v;
  UVector3 ps, delta, qa, qb, qc, qd;

  int face = -1;
  double normSign = outgoing ? 1 : -1;

  do
  {
    if (face == segment) continue;
    UPolyhedraSideVec& vec = (face == -1) ? vecs[segment] : vecs[face];
    //
    // Correct normal?
    //
    double dotProd = normSign * v.Dot(vec.normal);
    if (dotProd <= 0) continue;

    //
    // Is this face in front of the point along the trajectory?
    //
    delta = p - vec.center;
    distFromSurface = -normSign * delta.Dot(vec.normal);

    if (distFromSurface < -surfTolerance) continue;

    //
    //                            phi
    //      c -------- d           ^
    //      |          |           |
    //      a -------- b           +---> r/z
    //
    //
    // Do we remain on this particular segment?
    //
    qc = q - vec.edges[1]->corner[0];
    qd = q - vec.edges[1]->corner[1];

    if (normSign * qc.Cross(qd).Dot(v) < 0) continue;

    qa = q - vec.edges[0]->corner[0];
    qb = q - vec.edges[0]->corner[1];

    if (normSign * qa.Cross(qb).Dot(v) > 0) continue;

    //
    // We found the one and only segment we might be intersecting.
    // Do we remain within r/z bounds?
    //

    if (r[0] > 1 / UUtils::kInfinity && normSign * qa.Cross(qc).Dot(v) < 0) return false;
    if (r[1] > 1 / UUtils::kInfinity && normSign * qb.Cross(qd).Dot(v) > 0) return false;

    //
    // We allow the face to be slightly behind the trajectory
    // (surface tolerance) only if the point p is within
    // the vicinity of the face
    //
    if (distFromSurface < 0)
    {
      ps = p - vec.center;

      double rz = ps.Dot(vec.surfRZ);
      if (std::fabs(rz) > lenRZ + surfTolerance) return false;

      double pp = ps.Dot(vec.surfPhi);
      if (std::fabs(pp) > lenPhi[0] + lenPhi[1]*rz + surfTolerance) return false;
    }

    //
    // Intersection found. Return answer.
    //
    distance = distFromSurface / dotProd;
    normal = vec.normal;
    isAllBehind = allBehind;
    return true;
  }
  while (++face < numSide);

  //
  // Oh well. Better luck next time.
  //
  return false;
}


double UPolyhedraSide::Safety(const UVector3& p, bool outgoing)
{
  double normSign = outgoing ? -1 : +1;

  //
  // Try the closest phi segment first
  //
  int iPhi = ClosestPhiSegment(GetPhi(p));

  UVector3 pdotc = p - vecs[iPhi].center;
  double normDist = pdotc.Dot(vecs[iPhi].normal);

  if (normSign * normDist > -0.5 * VUSolid::Tolerance())
  {
    return DistanceAway(p, vecs[iPhi], &normDist);
  }

  //
  // Now we have an interesting problem... do we try to find the
  // closest facing side??
  //
  // Considered carefully, the answer is no. We know that if we
  // are asking for the distance out, we are supposed to be inside,
  // and vice versa.
  //

  return UUtils::kInfinity;
}


//
// Inside
//
VUSolid::EnumInside UPolyhedraSide::Inside(const UVector3& p,
                                           double tolerance,
                                           double* bestDistance)
{
  //
  // Which phi segment is closest to this point?
  //
  int iPhi = ClosestPhiSegment(GetPhi(p));

  double norm;
  //
  // Get distance to this segment
  //
  *bestDistance = DistanceToOneSide(p, vecs[iPhi], &norm);

  //
  // Use distance along normal to decide return value
  //
  if ((std::fabs(norm) < tolerance) && (*bestDistance < 2.0 * tolerance))
    return VUSolid::eSurface;

  if (norm < 0) return VUSolid::eInside;

  return VUSolid::eOutside;
}


//
// Normal
//
UVector3 UPolyhedraSide::Normal(const UVector3& p,
                                double* bestDistance)
{
  //
  // Which phi segment is closest to this point?
  //
  int iPhi = ClosestPhiSegment(GetPhi(p));

  //
  // Get distance to this segment
  //
  double norm;
  *bestDistance = DistanceToOneSide(p, vecs[iPhi], &norm);

  return vecs[iPhi].normal;
}


//
// Extent
//
double UPolyhedraSide::Extent(const UVector3 axis)
{
  if (axis.Perp2() < DBL_MIN)
  {
    //
    // Special case
    //
    return axis.z < 0 ? -cone->ZLo() : cone->ZHi();
  }

  int iPhi, i1, i2;
  double best;
  UVector3* list[4];

  //
  // Which phi segment, if any, does the axis belong to
  //
  iPhi = PhiSegment(GetPhi(axis));

  if (iPhi < 0)
  {
    //
    // No phi segment? Check front edge of first side and
    // last edge of second side
    //
    i1 = 0;
    i2 = numSide - 1;
  }
  else
  {
    //
    // Check all corners of matching phi side
    //
    i1 = iPhi;
    i2 = iPhi;
  }

  list[0] = vecs[i1].edges[0]->corner;
  list[1] = vecs[i1].edges[0]->corner + 1;
  list[2] = vecs[i2].edges[1]->corner;
  list[3] = vecs[i2].edges[1]->corner + 1;

  //
  // Who's biggest?
  //
  best = -UUtils::kInfinity;
  UVector3** vec = list;
  do
  {
    double answer = (*vec)->Dot(axis);
    if (answer > best) best = answer;
  }
  while (++vec < list + 4);

  return best;
}



//
// IntersectSidePlane
//
// Decide if a line correctly intersects one side plane of our segment.
// It is assumed that the correct side has been chosen, and thus only
// the z bounds (of the entire segment) are checked.
//
// normSign - To be multiplied against normal:
//            = +1.0 normal is unchanged
//            = -1.0 normal is reversed (now points inward)
//
// Arguments:
//  p    - (in) Point
//  v    - (in) Direction
//  vec    - (in) Description record of the side plane
//  normSign  - (in) Sign (+/- 1) to apply to normal
//  surfTolerance  - (in) Surface tolerance (generally > 0, see below)
//  distance  - (out) Distance along v to intersection
//  distFromSurface - (out) Distance from surface normal
//
// Notes:
//   surfTolerance  - Used to decide if a point is behind the surface,
//        a point is allow to be -surfTolerance behind the
//        surface (as measured along the normal), but *only*
//        if the point is within the r/z bounds + surfTolerance
//        of the segment.
//
bool UPolyhedraSide::IntersectSidePlane(const UVector3& p,
                                        const UVector3& v,
                                        const UPolyhedraSideVec& vec,
                                        double normSign,
                                        double surfTolerance,
                                        double& distance,
                                        double& distFromSurface)
{
  //
  // Correct normal? Here we have straight sides, and can safely ignore
  // intersections where the Dot product with the normal is zero.
  //
  double dotProd = normSign * v.Dot(vec.normal);

  if (dotProd <= 0) return false;

  //
  // Calculate distance to surface. If the side is too far
  // behind the point, we must reject it.
  //
  UVector3 delta = p - vec.center;
  distFromSurface = -normSign * delta.Dot(vec.normal);

  if (distFromSurface < -surfTolerance) return false;

  //
  // Calculate precise distance to intersection with the side
  // (along the trajectory, not normal to the surface)
  //
  distance = distFromSurface / dotProd;

  //
  // Do we fall off the r/z extent of the segment?
  //
  // Calculate this very, very carefully! Why?
  //         1. If a RZ end is at R=0, you can't miss!
  //         2. If you just fall off in RZ, the answer must
  //            be consistent with adjacent UPolyhedraSide faces.
  // (2) implies that only variables used by other UPolyhedraSide
  // faces may be used, which includes only: p, v, and the edge corners.
  // It also means that one side is a ">" or "<", which the other
  // must be ">=" or "<=". Fortunately, this isn't a new problem.
  // The solution below I borrowed from Joseph O'Rourke,
  // "Computational Geometry in C (Second Edition)"
  // See: http://cs.smith.edu/~orourke/
  //
  UVector3 ic = p + distance * v - vec.center;
  double atRZ = vec.surfRZ.Dot(ic);

  if (atRZ < 0)
  {
    if (r[0] == 0) return true;  // Can't miss!

    if (atRZ < -lenRZ * 1.2) return false; // Forget it! Missed by a mile.

    UVector3 q = p + v;
    UVector3 qa = q - vec.edges[0]->corner[0],
             qb = q - vec.edges[1]->corner[0];
    UVector3 qacb = qa.Cross(qb);
    if (normSign * qacb.Dot(v) < 0) return false;

    if (distFromSurface < 0)
    {
      if (atRZ < -lenRZ - surfTolerance) return false;
    }
  }
  else if (atRZ > 0)
  {
    if (r[1] == 0) return true;  // Can't miss!

    if (atRZ > lenRZ * 1.2) return false; // Missed by a mile

    UVector3 q = p + v;
    UVector3 qa = q - vec.edges[0]->corner[1],
             qb = q - vec.edges[1]->corner[1];
    UVector3 qacb = qa.Cross(qb);
    if (normSign * qacb.Dot(v) >= 0) return false;

    if (distFromSurface < 0)
    {
      if (atRZ > lenRZ + surfTolerance) return false;
    }
  }

  return true;
}


//
// LineHitsSegments
//
// Calculate which phi segments a line intersects in three dimensions.
// No check is made as to whether the intersections are within the z bounds of
// the segment.
//
int UPolyhedraSide::LineHitsSegments(const UVector3& p,
                                     const UVector3& v,
                                     int* i1, int* i2)
{
  double s1, s2;
  //
  // First, decide if and where the line intersects the cone
  //
  int n = cone->LineHitsCone(p, v, s1, s2);

  if (n == 0) return 0;

  //
  // Try first intersection.
  //
  *i1 = PhiSegment(std::atan2(p.y + s1 * v.y, p.x + s1 * v.x));
  if (n == 1)
  {
    return (*i1 < 0) ? 0 : 1;
  }

  //
  // Try second intersection
  //
  *i2 = PhiSegment(std::atan2(p.y + s2 * v.y, p.x + s2 * v.x));
  if (*i1 == *i2) return 0;

  if (*i1 < 0)
  {
    if (*i2 < 0) return 0;
    *i1 = *i2;
    return 1;
  }

  if (*i2 < 0) return 1;

  return 2;
}


//
// ClosestPhiSegment
//
// Decide which phi segment is closest in phi to the point.
// The result is the same as PhiSegment if there is no phi opening.
//
int UPolyhedraSide::ClosestPhiSegment(double phi0)
{
  int iPhi = PhiSegment(phi0);
  if (iPhi >= 0) return iPhi;

  //
  // Boogers! The points falls inside the phi segment.
  // Look for the closest point: the start, or  end
  //
  double phi = phi0;

  while (phi < startPhi) phi += 2 * UUtils::kPi;
  double d1 = phi - endPhi;

  while (phi > startPhi) phi -= 2 * UUtils::kPi;
  double d2 = startPhi - phi;

  return (d2 < d1) ? 0 : numSide - 1;
}


//
// PhiSegment
//
// Decide which phi segment an angle belongs to, counting from zero.
// A value of -1 indicates that the phi value is outside the shape
// (only possible if phiTotal < 360 degrees).
//
int UPolyhedraSide::PhiSegment(double phi0)
{
  //
  // How far are we from phiStart? Come up with a positive answer
  // that is less than 2*PI
  //
  double phi = phi0 - startPhi;
  while (phi < 0) phi += 2 * UUtils::kPi;
  while (phi > 2 * UUtils::kPi) phi -= 2 * UUtils::kPi;

  //
  // Divide
  //
  int answer = (int)(phi / deltaPhi);

  if (answer >= numSide)
  {
    if (phiIsOpen)
      return -1;  // Looks like we missed
    else
      answer = numSide - 1; // Probably just roundoff
  }

  return answer;
}


//
// GetPhi
//
// Calculate Phi for a given 3-vector (point), if not already cached for the
// same point, in the attempt to avoid consecutive computation of the same
// quantity
//
double UPolyhedraSide::GetPhi(const UVector3& p)
{
  double val = 0.;

  if (fPhi.first != p)
  {
    val = p.Phi();
    fPhi.first = p;
    fPhi.second = val;
  }
  else
  {
    val = fPhi.second;
  }
  return val;
}


//
// DistanceToOneSide
//
// Arguments:
//  p   - (in) Point to check
//  vec   - (in) vector Set of this side
//  normDist - (out) distance normal to the side or edge, as appropriate, signed
// Return value = total distance from the side
//
double UPolyhedraSide::DistanceToOneSide(const UVector3& p,
                                         const UPolyhedraSideVec& vec,
                                         double* normDist)
{
  UVector3 pct = p - vec.center;

  //
  // Get normal distance
  //
  *normDist = vec.normal.Dot(pct);

  //
  // Add edge penalty
  //
  return DistanceAway(p, vec, normDist);
}


//
// DistanceAway
//
// Add distance from side edges, if necesssary, to total distance,
// and updates normDist appropriate depending on edge normals.
//
double UPolyhedraSide::DistanceAway(const UVector3& p,
                                    const UPolyhedraSideVec& vec,
                                    double* normDist)
{
  double distOut2;
  UVector3 pct = p - vec.center;
  double distFaceNorm = *normDist;
  //
  // Okay, are we inside bounds?
  //
  double pcDotRZ  = pct.Dot(vec.surfRZ);
  double pcDotPhi = pct.Dot(vec.surfPhi);

  //
  // Go through all permutations.
  //                                                   Phi
  //               |              |                     ^
  //           B   |      H       |   E                 |
  //        ------[1]------------[3]-----               |
  //               |XXXXXXXXXXXXXX|                     +----> RZ
  //           C   |XXXXXXXXXXXXXX|   F
  //               |XXXXXXXXXXXXXX|
  //        ------[0]------------[2]----
  //           A   |      G       |   D
  //               |              |
  //
  // It's real messy, but at least it's quick
  //

  if (pcDotRZ < -lenRZ)
  {
    double lenPhiZ = lenPhi[0] - lenRZ * lenPhi[1];
    double distOutZ = pcDotRZ + lenRZ;
    //
    // Below in RZ
    //
    if (pcDotPhi < -lenPhiZ)
    {
      //
      // ...and below in phi. Find distance to point (A)
      //
      double distOutPhi = pcDotPhi + lenPhiZ;
      distOut2 = distOutPhi * distOutPhi + distOutZ * distOutZ;
      UVector3 pa = p - vec.edges[0]->corner[0];
      *normDist = pa.Dot(vec.edges[0]->cornNorm[0]);
    }
    else if (pcDotPhi > lenPhiZ)
    {
      //
      // ...and above in phi. Find distance to point (B)
      //
      double distOutPhi = pcDotPhi - lenPhiZ;
      distOut2 = distOutPhi * distOutPhi + distOutZ * distOutZ;
      UVector3 pb = p - vec.edges[1]->corner[0];
      *normDist = pb.Dot(vec.edges[1]->cornNorm[0]);
    }
    else
    {
      //
      // ...and inside in phi. Find distance to line (C)
      //
      UVector3 pa = p - vec.edges[0]->corner[0];
      distOut2 = distOutZ * distOutZ;
      *normDist = pa.Dot(vec.edgeNorm[0]);
    }
  }
  else if (pcDotRZ > lenRZ)
  {
    double lenPhiZ = lenPhi[0] + lenRZ * lenPhi[1];
    double distOutZ = pcDotRZ - lenRZ;
    //
    // Above in RZ
    //
    if (pcDotPhi < -lenPhiZ)
    {
      //
      // ...and below in phi. Find distance to point (D)
      //
      double distOutPhi = pcDotPhi + lenPhiZ;
      distOut2 = distOutPhi * distOutPhi + distOutZ * distOutZ;
      UVector3 pd = p - vec.edges[0]->corner[1];
      *normDist = pd.Dot(vec.edges[0]->cornNorm[1]);
    }
    else if (pcDotPhi > lenPhiZ)
    {
      //
      // ...and above in phi. Find distance to point (E)
      //
      double distOutPhi = pcDotPhi - lenPhiZ;
      distOut2 = distOutPhi * distOutPhi + distOutZ * distOutZ;
      UVector3 pe = p - vec.edges[1]->corner[1];
      *normDist = pe.Dot(vec.edges[1]->cornNorm[1]);
    }
    else
    {
      //
      // ...and inside in phi. Find distance to line (F)
      //
      distOut2 = distOutZ * distOutZ;
      UVector3 pd = p - vec.edges[0]->corner[1];
      *normDist = pd.Dot(vec.edgeNorm[1]);
    }
  }
  else
  {
    double lenPhiZ = lenPhi[0] + pcDotRZ * lenPhi[1];
    //
    // We are inside RZ bounds
    //
    if (pcDotPhi < -lenPhiZ)
    {
      //
      // ...and below in phi. Find distance to line (G)
      //
      double distOut = edgeNorm * (pcDotPhi + lenPhiZ);
      distOut2 = distOut * distOut;
      UVector3 pd = p - vec.edges[0]->corner[1];
      *normDist = pd.Dot(vec.edges[0]->normal);
    }
    else if (pcDotPhi > lenPhiZ)
    {
      //
      // ...and above in phi. Find distance to line (H)
      //
      double distOut = edgeNorm * (pcDotPhi - lenPhiZ);
      distOut2 = distOut * distOut;
      UVector3 pe = p - vec.edges[1]->corner[1];
      *normDist = pe.Dot(vec.edges[1]->normal);
    }
    else
    {
      //
      // Inside bounds! No penalty.
      //
      return std::fabs(distFaceNorm);
    }
  }
  return std::sqrt(distFaceNorm * distFaceNorm + distOut2);
}


//
// Calculation of surface area of a triangle.
// At the same time a random point in the triangle is given
//
double UPolyhedraSide::SurfaceTriangle(UVector3 p1,
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
// GetPointOnPlane
//
// Auxiliary method for GetPointOnSurface()
//
UVector3
UPolyhedraSide::GetPointOnPlane(UVector3 p0, UVector3 p1,
                                UVector3 p2, UVector3 p3,
                                double* Area)
{
  double chose, aOne, aTwo;
  UVector3 point1, point2;

  aOne = SurfaceTriangle(p0, p1, p2, &point1);
  aTwo = SurfaceTriangle(p2, p3, p0, &point2);
  *Area = aOne + aTwo;

  chose = UUtils::Random() * (aOne + aTwo);
  if ((chose >= 0.) && (chose < aOne))
  {
    return (point1);
  }
  return (point2);
}


//
// SurfaceArea()
//
double UPolyhedraSide::SurfaceArea()
{
  if (fSurfaceArea == 0.)
  {
    // Define the variables
    //
    double area, areas;
    UVector3 point1;
    UVector3 v1, v2, v3, v4;
    UPolyhedraSideVec* vec = vecs;
    areas = 0.;

    // Do a loop on all SideEdge
    //
    do
    {
      // Define 4points for a Plane or Triangle
      //
      v1 = vec->edges[0]->corner[0];
      v2 = vec->edges[0]->corner[1];
      v3 = vec->edges[1]->corner[1];
      v4 = vec->edges[1]->corner[0];
      point1 = GetPointOnPlane(v1, v2, v3, v4, &area);
      areas += area;
    }
    while (++vec < vecs + numSide);

    fSurfaceArea = areas;
  }
  return fSurfaceArea;
}


//
// GetPointOnFace()
//
UVector3 UPolyhedraSide::GetPointOnFace()
{
  // Define the variables
  //
  std::vector<double> areas;
  std::vector<UVector3> points;
  double area = 0;
  double result1;
  UVector3 point1;
  UVector3 v1, v2, v3, v4;
  UPolyhedraSideVec* vec = vecs;

  // Do a loop on all SideEdge
  //
  do
  {
    // Define 4points for a Plane or Triangle
    //
    v1 = vec->edges[0]->corner[0];
    v2 = vec->edges[0]->corner[1];
    v3 = vec->edges[1]->corner[1];
    v4 = vec->edges[1]->corner[0];
    point1 = GetPointOnPlane(v1, v2, v3, v4, &result1);
    points.push_back(point1);
    areas.push_back(result1);
    area += result1;
  }
  while (++vec < vecs + numSide);

  // Choose randomly one of the surfaces and point on it
  //
  double chose = area * UUtils::Random();
  double Achose1, Achose2;
  Achose1 = 0;
  Achose2 = 0.;
  int i = 0;
  do
  {
    Achose2 += areas[i];
    if (chose >= Achose1 && chose < Achose2)
    {
      point1 = points[i] ;
      break;
    }
    i++;
    Achose1 = Achose2;
  }
  while (i < numSide);

  return point1;
}
