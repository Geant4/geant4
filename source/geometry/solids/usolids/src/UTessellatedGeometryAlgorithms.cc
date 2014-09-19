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
// UTessellatedGeometryAlgorithms
//
// 11.07.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UTessellatedGeometryAlgorithms.hh"

#include <cfloat>

///////////////////////////////////////////////////////////////////////////////
//
// IntersectLineAndTriangle2D
//
// Determines whether there is an intersection between a line defined
// by r = p + s.v and a triangle defined by verticies p0, p0+e0 and p0+e1.
//
// Here:
//        p = 2D vector
//        s = scaler on [0,infinity)
//        v = 2D vector
//        p0, e0 and e1 are 2D vectors
// Information about where the intersection occurs is returned in the
// variable location.
//
// This is based on the work of Rickard Holmberg.
//
bool UTessellatedGeometryAlgorithms::IntersectLineAndTriangle2D(
  const UVector2& p,  const UVector2& v,
  const UVector2& p0, const UVector2& e0, const UVector2& e1,
  UVector2 location[2])
{
  UVector2 loc0[2];
  int e0i = IntersectLineAndLineSegment2D(p, v, p0, e0, loc0);
  if (e0i == 2)
  {
    location[0] = loc0[0];
    location[1] = loc0[1];
    return true;
  }

  UVector2 loc1[2];
  int e1i = IntersectLineAndLineSegment2D(p, v, p0, e1, loc1);
  if (e1i == 2)
  {
    location[0] = loc1[0];
    location[1] = loc1[1];
    return true;
  }

  if ((e0i == 1) && (e1i == 1))
  {
    if ((loc0[0] - p).mag2() < (loc1[0] - p).mag2())
    {
      location[0] = loc0[0];
      location[1] = loc1[0];
    }
    else
    {
      location[0] = loc1[0];
      location[1] = loc0[0];
    }
    return true;
  }

  UVector2 p1 = p0 + e0;
  UVector2 DE = e1 - e0;
  UVector2 loc2[2];
  int e2i = IntersectLineAndLineSegment2D(p, v, p1, DE, loc2);
  if (e2i == 2)
  {
    location[0] = loc2[0];
    location[1] = loc2[1];
    return true;
  }

  if ((e0i == 0) && (e1i == 0) && (e2i == 0)) return false;

  if ((e0i == 1) && (e2i == 1))
  {
    if ((loc0[0] - p).mag2() < (loc2[0] - p).mag2())
    {
      location[0] = loc0[0];
      location[1] = loc2[0];
    }
    else
    {
      location[0] = loc2[0];
      location[1] = loc0[0];
    }
    return true;
  }

  if ((e1i == 1) && (e2i == 1))
  {
    if ((loc1[0] - p).mag2() < (loc2[0] - p).mag2())
    {
      location[0] = loc1[0];
      location[1] = loc2[0];
    }
    else
    {
      location[0] = loc2[0];
      location[1] = loc1[0];
    }
    return true;
  }

  return false;
}

///////////////////////////////////////////////////////////////////////////////
//
// IntersectLineAndLineSegment2D
//
// Determines whether there is an intersection between a line defined
// by r = p0 + s.d0 and a line-segment with endpoints p1 and p1+d1.
// Here:
//        p0 = 2D vector
//        s  = scaler on [0,infinity)
//        d0 = 2D vector
//        p1 and d1 are 2D vectors
//
// This function returns:
// 0 - if there is no intersection;
// 1 - if there is a unique intersection;
// 2 - if the line and line-segments overlap, and the intersection is a
//     segment itself.
// Information about where the intersection occurs is returned in the
// as ??.
//
// This is based on the work of Rickard Holmberg as well as material published
// by Philip J Schneider and David H Eberly, "Geometric Tools for Computer
// Graphics," ISBN 1-55860-694-0, pp 244-245, 2003.
//
int UTessellatedGeometryAlgorithms::IntersectLineAndLineSegment2D(
  const UVector2& p0, const UVector2& d0,
  const UVector2& p1, const UVector2& d1,
  UVector2 location[2])
{
  UVector2 e     = p1 - p0;
  double kross    = Cross(d0, d1);
  double sqrKross = kross * kross;
  double sqrLen0  = d0.mag2();
  double sqrLen1  = d1.mag2();
  location[0]       = UVector2(0.0, 0.0);
  location[1]       = UVector2(0.0, 0.0);

  if (sqrKross > DBL_EPSILON * DBL_EPSILON * sqrLen0 * sqrLen1)
  {
    //
    //
    // The line and line segment are not parallel.  Determine if the intersection
    // is in positive s where r = p0 + s*d0, and for 0<=t<=1 where r = p1 + t*d1.
    //
    double ss = Cross(e, d1) / kross;
    if (ss < 0)          return 0; // Intersection does not occur for positive ss.
    double t = Cross(e, d0) / kross;
    if (t < 0 || t > 1) return 0; // Intersection does not occur on line-segment.
    //
    //
    // Intersection of lines is a single point on the forward-propagating line
    // defined by r = p0 + ss*d0, and the line segment defined by  r = p1 + t*d1.
    //
    location[0] = p0 + ss * d0;
    return 1;
  }
  //
  //
  // Line and line segment are parallel.  Determine whether they overlap or not.
  //
  double sqrLenE = e.mag2();
  kross            = Cross(e, d0);
  sqrKross         = kross * kross;
  if (sqrKross > DBL_EPSILON * DBL_EPSILON * sqrLen0 * sqrLenE)
  {
    return 0; //Lines are different.
  }
  //
  //
  // Lines are the same.  Test for overlap.
  //
  double s0   = d0.dot(e) / sqrLen0;
  double s1   = s0 + d0.dot(d1) / sqrLen0;
  double smin = 0.0;
  double smax = 0.0;

  if (s0 < s1)
  {
    smin = s0;
    smax = s1;
  }
  else
  {
    smin = s1;
    smax = s0;
  }

  if (smax < 0.0) return 0;
  else if (smin < 0.0)
  {
    location[0] = p0;
    location[1] = p0 + smax * d0;
    return 2;
  }
  else
  {
    location[0] = p0 + smin * d0;
    location[1] = p0 + smax * d0;
    return 2;
  }
}

double UTessellatedGeometryAlgorithms::Cross(const UVector2& v1,
                                             const UVector2& v2)
{
  return v1.x * v2.y - v1.y * v2.x;
}
