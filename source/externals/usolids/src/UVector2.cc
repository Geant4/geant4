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
// UVector2
//
// 19.09.12 Marek Gayer
//          Created from original implementation in CLHEP
// --------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include "UVector2.hh"

double UVector2::tolerance = UVector2::ZMpvToleranceTicks * 2.22045e-16;

double UVector2::setTolerance(double tol)
{
  // Set the tolerance for UVector2s to be considered near one another
  double oldTolerance(tolerance);
  tolerance = tol;
  return oldTolerance;
}

double UVector2::operator()(int i) const
{
  if (i == 0)
  {
    return x;
  }
  else if (i == 1)
  {
    return y;
  }
  else
  {
    // ZMthrowA(ZMxpvIndexRange("UVector2::operator(): bad index"));
    return 0.0;
  }
}

double& UVector2::operator()(int i)
{
  static double dummy;

  switch (i)
  {
    case X:
      return x;
    case Y:
      return y;
    default:
      // ZMthrowA (ZMxpvIndexRange("UVector2::operator() : bad index"));
      return dummy;
  }
}

void UVector2::rotate(double angler)
{
  double s = std::sin(angler);
  double c = std::cos(angler);
  double xx = x;
  x = c * xx - s * y;
  y = s * xx + c * y;
}

UVector2 operator/ (const UVector2& p, double a)
{
  if (a == 0)
  {
    // ZMthrowA(ZMxpvInfiniteVector( "Division of UVector2 by zero"));
  }
  return UVector2(p.x / a, p.y / a);
}

std::ostream& operator << (std::ostream& os, const UVector2& q)
{
  os << "(" << q.x << ", " << q.y << ")";
  return os;
}

//void ZMinput2doubles ( std::istream & is, const char * type,
//  double & x, double & y );

/*
std::istream & operator>>(std::istream & is, UVector2 & p) {
  double x, y;
  ZMinput2doubles ( is, "UVector2", x, y );
  p.set(x, y);
  return  is;
}  // operator>>()
*/

UVector2::operator UVector3() const
{
  return UVector3(x, y, 0.0);
}

int UVector2::compare(const UVector2& v) const
{
  if (y > v.y)
  {
    return 1;
  }
  else if (y < v.y)
  {
    return -1;
  }
  else if (x > v.x)
  {
    return 1;
  }
  else if (x < v.x)
  {
    return -1;
  }
  else
  {
    return 0;
  }
} /* Compare */


bool UVector2::operator > (const UVector2& v) const
{
  return (compare(v)  > 0);
}
bool UVector2::operator < (const UVector2& v) const
{
  return (compare(v)  < 0);
}
bool UVector2::operator>= (const UVector2& v) const
{
  return (compare(v) >= 0);
}
bool UVector2::operator<= (const UVector2& v) const
{
  return (compare(v) <= 0);
}

bool UVector2::isNear(const UVector2& p, double epsilon) const
{
  double limit = dot(p) * epsilon * epsilon;
  return ((*this - p).mag2() <= limit);
} /* isNear() */

double UVector2::howNear(const UVector2& p) const
{
  double d   = (*this - p).mag2();
  double pdp = dot(p);
  if ((pdp > 0) && (d < pdp))
  {
    return std::sqrt(d / pdp);
  }
  else if ((pdp == 0) && (d == 0))
  {
    return 0;
  }
  else
  {
    return 1;
  }
} /* howNear */

double UVector2::howParallel(const UVector2& v) const
{
  // | V1 x V2 | / | V1 dot V2 |
  // Of course, the "cross product" is fictitious but the math is valid
  double v1v2 = std::fabs(dot(v));
  if (v1v2 == 0)
  {
    // Zero is parallel to no other vector except for zero.
    return ((mag2() == 0) && (v.mag2() == 0)) ? 0 : 1;
  }
  double abscross = std::fabs(x * v.y - y - v.x);
  if (abscross >= v1v2)
  {
    return 1;
  }
  else
  {
    return abscross / v1v2;
  }
} /* howParallel() */

bool UVector2::isParallel(const UVector2& v,
                          double epsilon) const
{
  // | V1 x V2 | <= epsilon * | V1 dot V2 |
  // Of course, the "cross product" is fictitious but the math is valid
  double v1v2 = std::fabs(dot(v));
  if (v1v2 == 0)
  {
    // Zero is parallel to no other vector except for zero.
    return ((mag2() == 0) && (v.mag2() == 0));
  }
  double abscross = std::fabs(x * v.y - y - v.x);
  return (abscross <= epsilon * v1v2);
} /* isParallel() */

double UVector2::howOrthogonal(const UVector2& v) const
{
  // | V1 dot V2 | / | V1 x V2 |
  // Of course, the "cross product" is fictitious but the math is valid
  double v1v2 = std::fabs(dot(v));
  if (v1v2 == 0)
  {
    return 0; // Even if one or both are 0, they are considered orthogonal
  }
  double abscross = std::fabs(x * v.y - y - v.x);
  if (v1v2 >= abscross)
  {
    return 1;
  }
  else
  {
    return v1v2 / abscross;
  }
} /* howOrthogonal() */

bool UVector2::isOrthogonal(const UVector2& v,
                            double epsilon) const
{
  // | V1 dot V2 | <= epsilon * | V1 x V2 |
  // Of course, the "cross product" is fictitious but the math is valid
  double v1v2 = std::fabs(dot(v));
  double abscross = std::fabs(x * v.y - y - v.x);
  return (v1v2 <= epsilon * abscross);
} /* isOrthogonal() */

