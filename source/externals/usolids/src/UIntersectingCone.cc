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
// UIntersectingCone
//
// 19.02.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include <string>
#include <cmath>
#include <sstream>
#include "UIntersectingCone.hh"
#include "VUSolid.hh"

const double UIntersectingCone::EpsilonQuad = 1.0 / 9.0E99;

//
// Constructor
//
UIntersectingCone::UIntersectingCone(const double r[2],
                                     const double z[2])
{
  static const double halfCarTolerance
    = 0.5 * VUSolid::Tolerance(); // UGeometryTolerance::GetInstance()->GetSurfaceTolerance();

  //
  // What type of cone are we?
  //
  type1 = (std::fabs(z[1] - z[0]) > std::fabs(r[1] - r[0]));

  if (type1)
  {
    B = (r[1] - r[0]) / (z[1] - z[0]); // tube like
    A = 0.5 * (r[1] + r[0] - B * (z[1] + z[0]));
  }
  else
  {
    B = (z[1] - z[0]) / (r[1] - r[0]); // disk like
    A = 0.5 * (z[1] + z[0] - B * (r[1] + r[0]));
  }
  //
  // Calculate extent
  //
  if (r[0] < r[1])
  {
    rLo = r[0] - halfCarTolerance;
    rHi = r[1] + halfCarTolerance;
  }
  else
  {
    rLo = r[1] - halfCarTolerance;
    rHi = r[0] + halfCarTolerance;
  }

  if (z[0] < z[1])
  {
    zLo = z[0] - halfCarTolerance;
    zHi = z[1] + halfCarTolerance;
  }
  else
  {
    zLo = z[1] - halfCarTolerance;
    zHi = z[0] + halfCarTolerance;
  }
}


/*
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
UIntersectingCone::UIntersectingCone( __void__& )
: zLo(0.), zHi(0.), rLo(0.), rHi(0.), type1(false), A(0.), B(0.)
{
}
*/


//
// Destructor
//
UIntersectingCone::~UIntersectingCone()
{
}


//
// HitOn
//
// Check r or z extent, as appropriate, to see if the point is possibly
// on the cone.
//
bool UIntersectingCone::HitOn(const double r,
                              const double z)
{
  //
  // Be careful! The inequalities cannot be "<=" and ">=" here without
  // punching a tiny hole in our shape!
  //
  if (type1)
  {
    if (z < zLo || z > zHi) return false;
  }
  else
  {
    if (r < rLo || r > rHi) return false;
  }

  return true;
}


//
// LineHitsCone
//
// Calculate the intersection of a line with our conical surface, ignoring
// any phi division
//
int UIntersectingCone::LineHitsCone(const UVector3& p, const UVector3& v, double& s1, double& s2)
{
  if (type1)
  {
    return LineHitsCone1(p, v, s1, s2);
  }
  else
  {
    return LineHitsCone2(p, v, s1, s2);
  }
}


//
// LineHitsCone1
//
// Calculate the intersections of a line with a conical surface. Only
// suitable if zPlane[0] != zPlane[1].
//
// Equation of a line:
//
//       x = x0 + s*tx      y = y0 + s*ty     z = z0 + s*tz
//
// Equation of a conical surface:
//
//       x**2 + y**2 = (A + B*z)**2
//
// Solution is quadratic:
//
//  a*s**2 + b*s + c = 0
//
// where:
//
//  a = x0**2 + y0**2 - (A + B*z0)**2
//
//  b = 2*( x0*tx + y0*ty - (A*B - B*B*z0)*tz)
//
//  c = tx**2 + ty**2 - (B*tz)**2
//
// Notice, that if a < 0, this indicates that the two solutions (assuming
// they exist) are in opposite cones (that is, given z0 = -A/B, one z < z0
// and the other z > z0). For our shapes, the invalid solution is one
// which produces A + Bz < 0, or the one where Bz is smallest (most negative).
// Since Bz = B*s*tz, if B*tz > 0, we want the largest s, otherwise,
// the smaller.
//
// If there are two solutions on one side of the cone, we want to make
// sure that they are on the "correct" side, that is A + B*z0 + s*B*tz >= 0.
//
// If a = 0, we have a linear problem: s = c/b, which again gives one solution.
// This should be rare.
//
// For b*b - 4*a*c = 0, we also have one solution, which is almost always
// a line just grazing the surface of a the cone, which we want to ignore.
// However, there are two other, very rare, possibilities:
// a line intersecting the z axis and either:
//       1. At the same angle std::atan(B) to just miss one side of the cone, or
//       2. Intersecting the cone apex (0,0,-A/B)
// We *don't* want to miss these! How do we identify them? Well, since
// this case is rare, we can at least swallow a little more CPU than we would
// normally be comfortable with. Intersection with the z axis means
// x0*ty - y0*tx = 0. Case (1) means a==0, and we've already dealt with that
// above. Case (2) means a < 0.
//
// Now: x0*tx + y0*ty = 0 in terms of roundoff error. We can write:
//             Delta = x0*tx + y0*ty
//             b = 2*( Delta - (A*B + B*B*z0)*tz )
// For:
//             b*b - 4*a*c = epsilon
// where epsilon is small, then:
//             Delta = epsilon/2/B
//

/*
int UIntersectingCone::Solution (const UVector3 &p, const UVector3 &v, double a, double b, double c, double &s1, double &s2)
{
  return 0 || 1 || 2;
}
*/

int UIntersectingCone::LineHitsCone1(const UVector3& p, const UVector3& v, double& s1, double& s2)
{
  double x0 = p.x(), y0 = p.y(), z0 = p.z();
  double tx = v.x(), ty = v.y(), tz = v.z();

  double a = tx * tx + ty * ty - UUtils::sqr(B * tz);
  double b = 2 * (x0 * tx + y0 * ty - (A * B + B * B * z0) * tz);
  double c = x0 * x0 + y0 * y0 - UUtils::sqr(A + B * z0);

  double radical = b * b - 4 * a * c;
  double radicalSqrt;

  double minRadical = 1E-6 * std::fabs(b);

  if (radical < -minRadical)
  {
    return 0;  // No solution
  }

  if (radical < minRadical)
  {
    //
    // The radical is roughly zero: check for special, very rare, cases
    //
    if (std::fabs(a) > EpsilonQuad)
    {
      if (B == 0.)
      {
        return 0;
      }
      if (std::fabs(x0 * ty - y0 * tx) < std::fabs(1E-6 / B))
      {
        s1 = -0.5 * b / a;
        return 1;
      }
      return 0;
    }
    radicalSqrt = radical; //TODO: check this case
  }
  else
  {
    radicalSqrt = std::sqrt(radical);
  }

  if (a > EpsilonQuad)
  {
    double sa, sb, q = -0.5 * (b + (b < 0 ? -radicalSqrt : +radicalSqrt));
    sa = q / a;
    sb = c / q;
    if (sa < sb)
    {
      s1 = sa;
      s2 = sb;
    }
    else
    {
      s1 = sb;
      s2 = sa;
    }
    if (A + B * (z0 + (s1)*tz) < 0)
    {
      return 0;
    }
    // if ((z0 + (s1)*tz - A)/B < 0)  { return 0; } // these lines are equivalent
    return 2;
  }
  else if (a < -EpsilonQuad)
  {
    double sa, sb, q = -0.5 * (b + (b < 0 ? -radicalSqrt : +radicalSqrt));
    sa = q / a;
    sb = c / q;
    s1 = (tz * B > 0) ^ (sa > sb) ? sb : sa;
    return 1;
  }
  else if (std::fabs(b) < EpsilonQuad)
  {
    return 0;
  }
  else
  {
    s1 = -c / b;
    if (A + B * (z0 + (s1)*tz) < 0)
    {
      return 0;
    }
    return 1;
  }
}




int UIntersectingCone::LineHitsCone1Optimized(const UVector3& p, const UVector3& v, double& s1, double& s2)
{
  double x0 = p.x(), y0 = p.y(), z0 = p.z();
  double tx = v.x(), ty = v.y(), tz = v.z();

  double a = tx * tx + ty * ty - UUtils::sqr(B * tz);
  double b = 2 * (x0 * tx + y0 * ty - (A * B + B * B * z0) * tz);
  double c = x0 * x0 + y0 * y0 - UUtils::sqr(A + B * z0);

  double radical = b * b - 4 * a * c;

  double minRadical = 1E-6 * std::fabs(b);

  if (radical < -minRadical)
  {
    return 0;  // No solution
  }

  if (std::fabs(a) > EpsilonQuad)
  {
    if (radical < minRadical)
    {
      //
      // The radical is roughly zero: check for special, very rare, cases
      //

      if (B == 0.)
      {
        return 0;
      }
      if (std::fabs(x0 * ty - y0 * tx) < std::fabs(1E-6 / B))
      {
        s1 = -0.5 * b / a;
        return 1;
      }
      return 0;
    }
    else
    {
      double radicalSqrt = std::sqrt(radical);

      if (a > 0)
      {
        double sa, sb, q = -0.5 * (b + (b < 0 ? -radicalSqrt : +radicalSqrt));
        sa = q / a;
        sb = c / q;
        if (sa < sb)
        {
          s1 = sa;
          s2 = sb;
        }
        else
        {
          s1 = sb;
          s2 = sa;
        }
        if (A + B * (z0 + (s1)*tz) < 0)
        {
          return 0;
        }
        // if ((z0 + (s1)*tz - A)/B < 0)  { return 0; } // these lines are equivalent
        return 2;
      }
      else
      {
        double sa, sb, q = -0.5 * (b + (b < 0 ? -radicalSqrt : +radicalSqrt));
        sa = q / a;
        sb = c / q;
        s1 = (tz * B > 0) ^ (sa > sb) ? sb : sa;
        return 1;
      }

    }
  }

  if (std::fabs(b) < EpsilonQuad)
  {
    return 0;
  }
  else
  {
    s1 = -c / b;
    if (A + B * (z0 + (s1)*tz) < 0)
    {
      return 0;
    }
    return 1;
  }
}


//
// LineHitsCone2
//
// See comments under LineHitsCone1. In this routine, case2, we have:
//
//   Z = A + B*R
//
// The solution is still quadratic:
//
//  a = tz**2 - B*B*(tx**2 + ty**2)
//
//  b = 2*( (z0-A)*tz - B*B*(x0*tx+y0*ty) )
//
//  c = ( (z0-A)**2 - B*B*(x0**2 + y0**2) )
//
// The rest is much the same, except some details.
//
// a > 0 now means we intersect only once in the correct hemisphere.
//
// a > 0 ? We only want solution which produces R > 0.
// since R = (z0+s*tz-A)/B, for tz/B > 0, this is the largest s
//          for tz/B < 0, this is the smallest s
// thus, same as in case 1 ( since sign(tz/B) = sign(tz*B) )
//
int UIntersectingCone::LineHitsCone2(const UVector3& p,
                                     const UVector3& v,
                                     double& s1, double& s2)
{
  double x0 = p.x(), y0 = p.y(), z0 = p.z();
  double tx = v.x(), ty = v.y(), tz = v.z();

  // Special case which might not be so rare: B = 0 (precisely)
  //
  if (B == 0)
  {
    if (std::fabs(tz) < EpsilonQuad)
    {
      return 0;
    }

    s1 = (A - z0) / tz;
    return 1;
  }

  double B2 = B * B;

  double a = tz * tz - B2 * (tx * tx + ty * ty);
  double b = 2 * ((z0 - A) * tz - B2 * (x0 * tx + y0 * ty));
  double c = UUtils::sqr(z0 - A) - B2 * (x0 * x0 + y0 * y0);

  double radical = b * b - 4 * a * c;

  if (radical < -1E-6 * std::fabs(b))
  {
    return 0;  // No solution
  }

  if (radical < 1E-6 * std::fabs(b))
  {
    //
    // The radical is roughly zero: check for special, very rare, cases
    //
    if (std::fabs(a) > EpsilonQuad)
    {
      if (std::fabs(x0 * ty - y0 * tx) < std::fabs(1E-6 / B))
      {
        s1 = -0.5 * b / a;
        return 1;
      }
      return 0;
    }
  }
  else
  {
    radical = std::sqrt(radical);
  }

  if (a < -EpsilonQuad)
  {
    double sa, sb, q = -0.5 * (b + (b < 0 ? -radical : +radical));
    sa = q / a;
    sb = c / q;
    if (sa < sb)
    {
      s1 = sa;
      s2 = sb;
    }
    else
    {
      s1 = sb;
      s2 = sa;
    }
    if ((z0 + (s1)*tz - A) / B < 0)
    {
      return 0;
    }
    return 2;
  }
  else if (a > EpsilonQuad)
  {
    double sa, sb, q = -0.5 * (b + (b < 0 ? -radical : +radical));
    sa = q / a;
    sb = c / q;
    s1 = (tz * B > 0) ^ (sa > sb) ? sb : sa;
    return 1;
  }
  else if (std::fabs(b) < EpsilonQuad)
  {
    return 0;
  }
  else
  {
    s1 = -c / b;
    if ((z0 + (s1)*tz - A) / B < 0)
    {
      return 0;
    }
    return 1;
  }
}

