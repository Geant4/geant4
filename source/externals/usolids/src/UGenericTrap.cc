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
// UGenericTrap
//
// 21.10.13 Tatiana Nikitina, CERN; Ivana Hrivnacova, IPN Orsay
//          Adapted from Root Arb8 implementation
// --------------------------------------------------------------------

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

#include "UGenericTrap.hh"
#include "UUtils.hh"

#include "UTessellatedSolid.hh"
#include "UTriangularFacet.hh"
#include "UQuadrangularFacet.hh"

const int    UGenericTrap::fgkNofVertices = 8;
const double UGenericTrap::fgkTolerance = 1E-3;

// --------------------------------------------------------------------

UGenericTrap::UGenericTrap(const std::string& name, double halfZ,
                           const std::vector<UVector2>&  vertices)
  : VUSolid(name),
    fDz(halfZ),
    fVertices(),
    fIsTwisted(false),
    fTessellatedSolid(0),
    fMinBBoxVector(UVector3(0, 0, 0)),
    fMaxBBoxVector(UVector3(0, 0, 0)),
    fVisSubdivisions(0),
    fBoundBox(0),
    fSurfaceArea(0.),
    fCubicVolume(0.)

{
  // General constructor
  Initialise(vertices);
}

// --------------------------------------------------------------------

UGenericTrap::UGenericTrap()
  : VUSolid(""),
    fDz(0.),
    fVertices(),
    fIsTwisted(false),
    fTessellatedSolid(0),
    fMinBBoxVector(UVector3(0,0,0)),
    fMaxBBoxVector(UVector3(0,0,0)),
    fVisSubdivisions(0),
    fBoundBox(0),
    fSurfaceArea(0.),
    fCubicVolume(0.)
{
  // Fake default constructor - sets only member data and allocates memory
  //                            for usage restricted to object persistency.
}

// --------------------------------------------------------------------

UGenericTrap::~UGenericTrap()
{
  // Destructor
  delete fTessellatedSolid;
  delete fBoundBox;
}

// --------------------------------------------------------------------

UGenericTrap::UGenericTrap(const UGenericTrap& rhs)
  : VUSolid(rhs),
    fDz(rhs.fDz), fVertices(rhs.fVertices),
    fIsTwisted(rhs.fIsTwisted), fTessellatedSolid(0),
    fMinBBoxVector(rhs.fMinBBoxVector), fMaxBBoxVector(rhs.fMaxBBoxVector),
    fVisSubdivisions(rhs.fVisSubdivisions), fBoundBox(0),
    fSurfaceArea(rhs.fSurfaceArea), fCubicVolume(rhs.fCubicVolume)
{
   for (size_t i=0; i<4; ++i)  { fTwist[i] = rhs.fTwist[i]; }
   ComputeBBox();
#ifdef UTESS_TEST
   if (rhs.fTessellatedSolid && !fIsTwisted )
   { fTessellatedSolid = CreateTessellatedSolid(); }
#endif
}

// --------------------------------------------------------------------

UGenericTrap& UGenericTrap::operator = (const UGenericTrap& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   VUSolid::operator=(rhs);

   // Copy data
   //
   fDz = rhs.fDz; fVertices = rhs.fVertices;
   fIsTwisted = rhs.fIsTwisted; fTessellatedSolid = 0;
   fMinBBoxVector = rhs.fMinBBoxVector; fMaxBBoxVector = rhs.fMaxBBoxVector;
   fVisSubdivisions = rhs.fVisSubdivisions;
   fSurfaceArea = rhs.fSurfaceArea; fCubicVolume = rhs.fCubicVolume;
   for (size_t i=0; i<4; ++i)  { fTwist[i] = rhs.fTwist[i]; }
   delete fBoundBox; ComputeBBox();
#ifdef UTESS_TEST
   if (rhs.fTessellatedSolid && !fIsTwisted )
   { delete fTessellatedSolid; fTessellatedSolid = CreateTessellatedSolid(); }
#endif

   return *this;
}

// --------------------------------------------------------------------

void UGenericTrap::Initialise(const std::vector<UVector2>&  vertices)
{
  const double min_length = 5 * 1.e-6;
  double  length = 0.;
  int  k = 0;
  std::string errorDescription = "InvalidSetup in \" ";
  errorDescription += GetName();
  errorDescription += "\"";

  // Check vertices size

  if (int (vertices.size()) != fgkNofVertices)
  {
    UUtils::Exception("UGenericTrap::UGenericTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, "Number of vertices != 8");
  }

  // Check dZ
  //
  if (fDz < VUSolid::fgTolerance)
  {
    UUtils::Exception("UGenericTrap::UGenericTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, "dZ is too small or negative");
  }

  // Check Ordering and Copy vertices
  //
  if (CheckOrder(vertices))
  {
    for (int  i = 0; i < fgkNofVertices; ++i)
    {
      fVertices.push_back(vertices[i]);
    }
  }
  else
  {
    for (int  i = 0; i < 4; ++i)
    {
      fVertices.push_back(vertices[3 - i]);
    }
    for (int  i = 0; i < 4; ++i)
    {
      fVertices.push_back(vertices[7 - i]);
    }
  }

  // Check length of segments and Adjust
  //
  for (int  j = 0; j < 2; j++)
  {
    for (int  i = 1; i < 4; ++i)
    {
      k = j * 4 + i;
      length = (fVertices[k] - fVertices[k - 1]).mag();
      if ((length < min_length) && (length > VUSolid::fgTolerance))
      {
        std::ostringstream message;
        message << "Length segment is too small." << std::endl
                << "Distance between " << fVertices[k - 1] << " and "
                << fVertices[k] << " is only " << length << " mm !"
                << "Vertices will be collapsed.";
        UUtils::Exception("UGenericTrap::UGenericTrap()", "GeomSolids1001",
                          UWarning, 1, message.str().c_str());
        fVertices[k] = fVertices[k - 1];
      }
    }
  }

  // Compute Twist
  //
  for (int  i = 0; i < 4; i++)
  {
    fTwist[i] = 0.;
  }
  fIsTwisted = ComputeIsTwisted();

  // Compute Bounding Box
  //
  ComputeBBox();

  // If not twisted - create tessellated solid
  // (an alternative implementation for testing)
  //
#ifdef UTESS_TEST
  if (!fIsTwisted)
  {
    fTessellatedSolid = CreateTessellatedSolid();
  }
#endif
}

// --------------------------------------------------------------------

VUSolid::EnumInside
UGenericTrap::InsidePolygone(const UVector3& p, const UVector2* poly) const
{
  static const double  halfCarTolerance = VUSolid::fgTolerance * 0.5;
  VUSolid::EnumInside  in = eInside;
  double  cross, len2;
  int  count = 0;
  UVector2 pt(p.x(), p.y());

  for (int  i = 0; i < 4; i++)
  {
    int  j = (i + 1) % 4;

    cross = (p.x() - poly[i].x) * (poly[j].y - poly[i].y) -
            (p.y() - poly[i].y) * (poly[j].x - poly[i].x);
    if (cross < 0.)return eOutside;

    len2 = (poly[i] - poly[j]).mag2();
    if (len2 > VUSolid::fgTolerance)
    {
      if (cross * cross <= len2 * halfCarTolerance * halfCarTolerance) // Surface check
      {
        double  test;
        int  k, l;
        if (poly[i].y > poly[j].y)
        {
          k = i;
          l = j;
        }
        else
        {
          k = j;
          l = i;
        }

        if (poly[k].x != poly[l].x)
        {
          test = (p.x() - poly[l].x) / (poly[k].x - poly[l].x)
                 * (poly[k].y - poly[l].y) + poly[l].y;
        }
        else
        {
          test = p.y();
        }

        // Check if point is Inside Segment
        //
        if ((test >= (poly[l].y - halfCarTolerance))
            && (test <= (poly[k].y + halfCarTolerance)))
        {
          return eSurface;
        }
        else
        {
          return eOutside;
        }
      }
      else if (cross < 0.)
      {
        return eOutside;
      }
    }
    else
    {
      count++;
    }
  }

  // All collapsed vertices, Tet like
  //
  if (count == 4)
  {
    if ((std::fabs(p.x() - poly[0].x) + std::fabs(p.y() - poly[0].y)) > halfCarTolerance)
    {
      in = eOutside;
    }
  }
  return in;


}

//
//_____________________________________________________________________________

bool UGenericTrap::IsSameLine(const UVector2& p,
                              const UVector2& l1, const UVector2& l2) const
{
  // Return true if p is on the line through l1, l2

  if (l1.x == l2.x)
  {
    return std::fabs(p.x - l1.x) < VUSolid::fgTolerance * 0.5;
  }
  double  slope = ((l2.y - l1.y) / (l2.x - l1.x));
  double predy = l1.y +  slope * (p.x - l1.x);
  double dy = p.y - predy;

  // Calculate perpendicular distance
  //
  // double perpD= std::fabs(dy) / std::sqrt( 1 + slope * slope );
  // bool   simpleComp= (perpD<0.5*VUSolid::fgTolerance);

  // Check perpendicular distance vs tolerance 'directly'
  //
  const double tol = 0.5 * VUSolid::fgTolerance ;
  bool    squareComp = (dy * dy < (1 + slope * slope) * tol * tol);

  // return  simpleComp;
  return squareComp;
}

//_____________________________________________________________________________

bool UGenericTrap::IsSameLineSegment(const UVector2& p,
                                     const UVector2& l1, const UVector2& l2) const
{
  // Return true if p is on the line through l1, l2 and lies between
  // l1 and l2

  if (p.x < std::min(l1.x, l2.x) - VUSolid::fgTolerance * 0.5 ||
      p.x > std::max(l1.x, l2.x) + VUSolid::fgTolerance * 0.5 ||
      p.y < std::min(l1.y, l2.y) - VUSolid::fgTolerance * 0.5 ||
      p.y > std::max(l1.y, l2.y) + VUSolid::fgTolerance * 0.5)
  {
    return false;
  }

  return IsSameLine(p, l1, l2);
}

// --------------------------------------------------------------------

VUSolid::EnumInside UGenericTrap::Inside(const UVector3& p) const
{
  // Test if point is inside this shape

#ifdef G4TESS_TEST
  if (fTessellatedSolid)
  {
    return fTessellatedSolid->Inside(p);
  }
#endif

  static const double  halfCarTolerance = VUSolid::fgTolerance * 0.5;
  VUSolid::EnumInside innew = eOutside;
  UVector2 xy[4];
  if (fBoundBox->Inside(p) == eOutside) return eOutside;

  if (std::fabs(p.z()) <= fDz + halfCarTolerance) // First check Z range
  {
    // Compute intersection between Z plane containing point and the shape
    //
    double  cf = 0.5 * (fDz - p.z()) / fDz;
    for (int  i = 0; i < 4; i++)
    {
      xy[i] = fVertices[i + 4] + cf * (fVertices[i] - fVertices[i + 4]);
    }

    innew = InsidePolygone(p, xy);

    if ((innew == eInside) || (innew == eSurface))
    {
      if (std::fabs(p.z()) > fDz - halfCarTolerance)
      {
        innew = eSurface;
      }
    }

  }
  return innew;
}
// --------------------------------------------------------------------

bool UGenericTrap::Normal(const UVector3& p, UVector3& aNormal) const
{
  // Calculate side nearest to p, and return normal
  // If two sides are equidistant, sum of the Normal is returned

#ifdef G4TESS_TEST
  if (fTessellatedSolid)
  {
    //return fTessellatedSolid->Normal(p);
  }
#endif

  static const double  halfCarTolerance = VUSolid::fgTolerance * 0.5;

  UVector3 lnorm, sumnorm(0., 0., 0.), apprnorm(0., 0., 1.),
           p0, p1, p2, r1, r2, r3, r4;
  int  noSurfaces = 0;
  double  distxy, distz;
  bool zPlusSide = false;

  distz = fDz - std::fabs(p.z());
  if (distz < halfCarTolerance)
  {
    if (p.z() > 0)
    {
      zPlusSide = true;
      sumnorm = UVector3(0, 0, 1);
    }
    else
    {
      sumnorm = UVector3(0, 0, -1);
    }
    noSurfaces ++;
  }

  // Check lateral planes
  //
  UVector2 vertices[4];
  double  cf = 0.5 * (fDz - p.z()) / fDz;
  for (int  i = 0; i < 4; i++)
  {
    vertices[i] = fVertices[i + 4] + cf * (fVertices[i] - fVertices[i + 4]);
  }


  // Compute distance for lateral planes
  //
  for (int  q = 0; q < 4; q++)
  {
    p0 = UVector3(vertices[q].x, vertices[q].y, p.z());
    if (zPlusSide)
    {
      p1 = UVector3(fVertices[q].x, fVertices[q].y, -fDz);
    }
    else
    {
      p1 = UVector3(fVertices[q + 4].x, fVertices[q + 4].y, fDz);
    }
    p2 = UVector3(vertices[(q + 1) % 4].x, vertices[(q + 1) % 4].y, p.z());

    // Collapsed vertices
    //
    if ((p2 - p0).Mag2() < VUSolid::fgTolerance)
    {
      if (std::fabs(p.z() + fDz) > VUSolid::fgTolerance)
      {
        p2 = UVector3(fVertices[(q + 1) % 4].x, fVertices[(q + 1) % 4].y, -fDz);
      }
      else
      {
        p2 = UVector3(fVertices[(q + 1) % 4 + 4].x, fVertices[(q + 1) % 4 + 4].y, fDz);
      }
    }
    lnorm = (p1 - p0).Cross(p2 - p0);
    lnorm = lnorm.Unit();
    if (zPlusSide)
    {
      lnorm = -lnorm;
    }

    // Adjust Normal for Twisted Surface
    //
    if ((fIsTwisted) && (GetTwistAngle(q) != 0))
    {
      double  normP = (p2 - p0).Mag();
      if (normP)
      {
        double  proj = (p - p0).Dot(p2 - p0) / normP;
        if (proj < 0)
        {
          proj = 0;
        }
        if (proj > normP)
        {
          proj = normP;
        }
        int  j = (q + 1) % 4;
        r1 = UVector3(fVertices[q + 4].x, fVertices[q + 4].y, fDz);
        r2 = UVector3(fVertices[j + 4].x, fVertices[j + 4].y, fDz);
        r3 = UVector3(fVertices[q].x, fVertices[q].y, -fDz);
        r4 = UVector3(fVertices[j].x, fVertices[j].y, -fDz);
        r1 = r1 + proj * (r2 - r1) / normP;
        r3 = r3 + proj * (r4 - r3) / normP;
        r2 = r1 - r3;
        r4 = r2.Cross(p2 - p0);
        r4 = r4.Unit();
        lnorm = r4;
      }
    }   // End if fIsTwisted

    distxy = std::fabs((p0 - p).Dot(lnorm));
    if (distxy < halfCarTolerance)
    {
      noSurfaces ++;

      // Negative sign for Normal is taken for Outside Normal
      //
      sumnorm = sumnorm + lnorm;
    }

    // For ApproxSurfaceNormal
    //
    if (distxy < distz)
    {
      distz = distxy;
      apprnorm = lnorm;
    }
  }  // End for loop

  // Calculate final Normal, add Normal in the Corners and Touching Sides
  //
  if (noSurfaces == 0)
  {
    #ifdef UDEBUG
    UUtils::Exception("UGenericTrap::SurfaceNormal(p)", "GeomSolids1002",
                      UWarning, 1, "Point p is not on surface !?");
    if (Inside(p) != eSurface)
    {
      std::cout << "Point is not on Surface, confirmed by Inside()" << p << std::endl;
    }
    #endif
    sumnorm = apprnorm;
    // Add Approximative Surface Normal Calculation
  }
  else if (noSurfaces == 1)
  {
    sumnorm = sumnorm;
  }
  else
  {
    sumnorm = sumnorm.Unit();
  }

  aNormal = sumnorm;
  return  noSurfaces != 0 ;
}

// --------------------------------------------------------------------

UVector3 UGenericTrap::NormalToPlane(const UVector3& p,
                                     const int  ipl) const
{
  // Return normal to given lateral plane ipl

#ifdef G4TESS_TEST
  if (fTessellatedSolid)
  {
    return fTessellatedSolid->SurfaceNormal(p);
  }
#endif

  static const double  halfCarTolerance = VUSolid::fgTolerance * 0.5;
  UVector3 lnorm, norm(0., 0., 0.), p0, p1, p2;

  double   distz = fDz - p.z();
  int  i = ipl; // current plane index

  UVector2 u, v;
  UVector3 r1, r2, r3, r4;
  double  cf = 0.5 * (fDz - p.z()) / fDz;
  int  j = (i + 1) % 4;

  u = fVertices[i + 4] + cf * (fVertices[i] - fVertices[i + 4]);
  v = fVertices[j + 4] + cf * (fVertices[j] - fVertices[j + 4]);

  // Compute cross product
  //
  p0 = UVector3(u.x, u.y, p.z());

  if (std::fabs(distz) < halfCarTolerance)
  {
    p1 = UVector3(fVertices[i].x, fVertices[i].y, -fDz);
    distz = -1;
  }
  else
  {
    p1 = UVector3(fVertices[i + 4].x, fVertices[i + 4].y, fDz);
  }
  p2 = UVector3(v.x, v.y, p.z());

  // Collapsed vertices
  //
  if ((p2 - p0).Mag2() < VUSolid::fgTolerance)
  {
    if (std::fabs(p.z() + fDz) > halfCarTolerance)
    {
      p2 = UVector3(fVertices[j].x, fVertices[j].y, -fDz);
    }
    else
    {
      p2 = UVector3(fVertices[j + 4].x, fVertices[j + 4].y, fDz);
    }
  }
  lnorm = -(p1 - p0).Cross(p2 - p0);
  if (distz > -halfCarTolerance)
  {
    lnorm = -lnorm.Unit();
  }
  else
  {
    lnorm = lnorm.Unit();
  }

  // Adjust Normal for Twisted Surface
  //
  if ((fIsTwisted) && (GetTwistAngle(ipl) != 0))
  {
    double  normP = (p2 - p0).Mag();
    if (normP)
    {
      double  proj = (p - p0).Dot(p2 - p0) / normP;
      if (proj < 0)
      {
        proj = 0;
      }
      if (proj > normP)
      {
        proj = normP;
      }

      r1 = UVector3(fVertices[i + 4].x, fVertices[i + 4].y, fDz);
      r2 = UVector3(fVertices[j + 4].x, fVertices[j + 4].y, fDz);
      r3 = UVector3(fVertices[i].x, fVertices[i].y, -fDz);
      r4 = UVector3(fVertices[j].x, fVertices[j].y, -fDz);
      r1 = r1 + proj * (r2 - r1) / normP;
      r3 = r3 + proj * (r4 - r3) / normP;
      r2 = r1 - r3;
      r4 = r2.Cross(p2 - p0);
      r4 = r4.Unit();
      lnorm = r4;
    }
  }  // End if fIsTwisted

  return lnorm;
}

// --------------------------------------------------------------------

double  UGenericTrap::DistToPlane(const UVector3& p,
                                  const UVector3& v,
                                  const int  ipl) const
{
  // Computes distance to plane ipl :
  // ipl=0 : points 0,4,1,5
  // ipl=1 : points 1,5,2,6
  // ipl=2 : points 2,6,3,7
  // ipl=3 : points 3,7,0,4

  static const double  halfCarTolerance = 0.5 * VUSolid::fgTolerance;
  double  xa, xb, xc, xd, ya, yb, yc, yd;

  int  j = (ipl + 1) % 4;

  xa = fVertices[ipl].x;
  ya = fVertices[ipl].y;
  xb = fVertices[ipl + 4].x;
  yb = fVertices[ipl + 4].y;
  xc = fVertices[j].x;
  yc = fVertices[j].y;
  xd = fVertices[4 + j].x;
  yd = fVertices[4 + j].y;

  double  dz2 = 0.5 / fDz;
  double  tx1 = dz2 * (xb - xa);
  double  ty1 = dz2 * (yb - ya);
  double  tx2 = dz2 * (xd - xc);
  double  ty2 = dz2 * (yd - yc);
  double  dzp = fDz + p.z();
  double  xs1 = xa + tx1 * dzp;
  double  ys1 = ya + ty1 * dzp;
  double  xs2 = xc + tx2 * dzp;
  double  ys2 = yc + ty2 * dzp;
  double  dxs = xs2 - xs1;
  double  dys = ys2 - ys1;
  double  dtx = tx2 - tx1;
  double  dty = ty2 - ty1;

  double  a = (dtx * v.y() - dty * v.x() + (tx1 * ty2 - tx2 * ty1) * v.z()) * v.z();
  double  b = dxs * v.y() - dys * v.x() + (dtx * p.y() - dty * p.x() + ty2 * xs1 - ty1 * xs2
                                       + tx1 * ys2 - tx2 * ys1) * v.z();
  double  c = dxs * p.y() - dys * p.x() + xs1 * ys2 - xs2 * ys1;
  double  q = UUtils::kInfinity;
  double  x1, x2, y1, y2, xp, yp, zi;

  if (std::fabs(a) < VUSolid::fgTolerance)
  {
    if (std::fabs(b) < VUSolid::fgTolerance)
    {
      return UUtils::kInfinity;
    }
    q = -c / b;

    // Check if Point is on the Surface

    if (q > -halfCarTolerance)
    {
      if (q < halfCarTolerance)
      {
        if (NormalToPlane(p, ipl).Dot(v) <= 0)
        {
          if (Inside(p) != eOutside)
          {
            return 0.;
          }
        }
        else
        {
          return UUtils::kInfinity;
        }
      }

      // Check the Intersection
      //
      zi = p.z() + q * v.z();
      if (std::fabs(zi) < fDz)
      {
        x1 = xs1 + tx1 * v.z() * q;
        x2 = xs2 + tx2 * v.z() * q;
        xp = p.x() + q * v.x();
        y1 = ys1 + ty1 * v.z() * q;
        y2 = ys2 + ty2 * v.z() * q;
        yp = p.y() + q * v.y();
        zi = (xp - x1) * (xp - x2) + (yp - y1) * (yp - y2);
        if (zi <= halfCarTolerance)
        {
          return q;
        }
      }
    }
    return UUtils::kInfinity;
  }
  double  d = b * b - 4 * a * c;
  if (d >= 0)
  {
    if (a > 0)
    {
      q = 0.5 * (-b - std::sqrt(d)) / a;
    }
    else
    {
      q = 0.5 * (-b + std::sqrt(d)) / a;
    }

    // Check if Point is on the Surface
    //
    if (q > -halfCarTolerance)
    {
      if (q < halfCarTolerance)
      {
        if (NormalToPlane(p, ipl).Dot(v) <= 0)
        {
          if (Inside(p) != eOutside)
          {
            return 0.;
          }
        }
        else  // Check second root; return UUtils::kInfinity
        {
          if (a > 0)
          {
            q = 0.5 * (-b + std::sqrt(d)) / a;
          }
          else
          {
            q = 0.5 * (-b - std::sqrt(d)) / a;
          }
          if (q <= halfCarTolerance)
          {
            return UUtils::kInfinity;
          }
        }
      }
      // Check the Intersection
      //
      zi = p.z() + q * v.z();
      if (std::fabs(zi) < fDz)
      {
        x1 = xs1 + tx1 * v.z() * q;
        x2 = xs2 + tx2 * v.z() * q;
        xp = p.x() + q * v.x();
        y1 = ys1 + ty1 * v.z() * q;
        y2 = ys2 + ty2 * v.z() * q;
        yp = p.y() + q * v.y();
        zi = (xp - x1) * (xp - x2) + (yp - y1) * (yp - y2);
        if (zi <= halfCarTolerance)
        {
          return q;
        }
      }
    }
    if (a > 0)
    {
      q = 0.5 * (-b + std::sqrt(d)) / a;
    }
    else
    {
      q = 0.5 * (-b - std::sqrt(d)) / a;
    }

    // Check if Point is on the Surface
    //
    if (q > -halfCarTolerance)
    {
      if (q < halfCarTolerance)
      {
        if (NormalToPlane(p, ipl).Dot(v) <= 0)
        {
          if (Inside(p) != eOutside)
          {
            return 0.;
          }
        }
        else   // Check second root; return UUtils::kInfinity.
        {
          if (a > 0)
          {
            q = 0.5 * (-b - std::sqrt(d)) / a;
          }
          else
          {
            q = 0.5 * (-b + std::sqrt(d)) / a;
          }
          if (q <= halfCarTolerance)
          {
            return UUtils::kInfinity;
          }
        }
      }
      // Check the Intersection
      //
      zi = p.z() + q * v.z();
      if (std::fabs(zi) < fDz)
      {
        x1 = xs1 + tx1 * v.z() * q;
        x2 = xs2 + tx2 * v.z() * q;
        xp = p.x() + q * v.x();
        y1 = ys1 + ty1 * v.z() * q;
        y2 = ys2 + ty2 * v.z() * q;
        yp = p.y() + q * v.y();
        zi = (xp - x1) * (xp - x2) + (yp - y1) * (yp - y2);
        if (zi <= halfCarTolerance)
        {
          return q;
        }
      }
    }
  }
  return UUtils::kInfinity;
}

// --------------------------------------------------------------------

double  UGenericTrap::DistanceToIn(const UVector3& p,
                                   const UVector3& v, double) const
{
#ifdef G4TESS_TEST
  if (fTessellatedSolid)
  {
    return fTessellatedSolid->DistanceToIn(p, v);
  }
#endif
  if (fBoundBox->DistanceToIn(p, v) == UUtils::kInfinity)return UUtils::kInfinity;
  static const double  halfCarTolerance = VUSolid::fgTolerance * 0.5;

  double  dist[5];
  UVector3 n;

  // Check lateral faces
  //
  int  i;
  for (i = 0; i < 4; i++)
  {
    dist[i] = DistToPlane(p, v, i);
  }

  // Check Z planes
  //
  dist[4] = UUtils::kInfinity;
  if (std::fabs(p.z()) > fDz - halfCarTolerance)
  {
    if (v.z())
    {
      UVector3 pt;
      if (p.z() > 0)
      {
        dist[4] = (fDz - p.z()) / v.z();
      }
      else
      {
        dist[4] = (-fDz - p.z()) / v.z();
      }
      if (dist[4] < -halfCarTolerance)
      {
        dist[4] = UUtils::kInfinity;
      }
      else
      {
        if (dist[4] < halfCarTolerance)
        {
          if (p.z() > 0)
          {
            n = UVector3(0, 0, 1);
          }
          else
          {
            n = UVector3(0, 0, -1);
          }
          if (n.Dot(v) < 0)
          {
            dist[4] = 0.;
          }
          else
          {
            dist[4] = UUtils::kInfinity;
          }
        }
        pt = p + dist[4] * v;
        if (Inside(pt) == eOutside)
        {
          dist[4] = UUtils::kInfinity;
        }
      }
    }
  }
  double  distmin = dist[0];
  for (i = 1; i < 5; i++)
  {
    if (dist[i] < distmin)
    {
      distmin = dist[i];
    }
  }

  if (distmin < halfCarTolerance)
  {
    distmin = 0.;
  }

  return distmin;
}

// --------------------------------------------------------------------

double UGenericTrap::SafetyFromOutside(const UVector3& p, bool precise) const
{
  // Computes the closest distance from given point to this shape

#ifdef G4TESS_TEST
  if (fTessellatedSolid)
  {
    return fTessellatedSolid->DistanceToIn(p);
  }
#endif
  if (!precise)return fBoundBox->SafetyFromOutside(p, !precise);
  double  safz = std::fabs(p.z()) - fDz;
  if (safz < 0)
  {
    safz = 0;
  }

  int  iseg;
  double  safe  = safz;
  double  safxy = safz;

  for (iseg = 0; iseg < 4; iseg++)
  {
    safxy = SafetyToFace(p, iseg);
    if (safxy > safe)
    {
      safe = safxy;
    }
  }

  return safe;
}

// --------------------------------------------------------------------

double
UGenericTrap::SafetyToFace(const UVector3& p, const int  iseg) const
{
  // Estimate distance to lateral plane defined by segment iseg in range [0,3]
  // Might be negative: plane seen only from inside

  UVector3 p1, norm;
  double  safe;

  p1 = UVector3(fVertices[iseg].x, fVertices[iseg].y, -fDz);

  norm = NormalToPlane(p, iseg);
  safe = (p - p1).Dot(norm); // Can be negative

  return safe;
}

// --------------------------------------------------------------------

double
UGenericTrap::DistToTriangle(const UVector3& p,
                             const UVector3& v, const int  ipl) const
{
  static const double  halfCarTolerance = VUSolid::fgTolerance * 0.5;

  double  xa = fVertices[ipl].x;
  double  ya = fVertices[ipl].y;
  double  xb = fVertices[ipl + 4].x;
  double  yb = fVertices[ipl + 4].y;
  int  j = (ipl + 1) % 4;
  double  xc = fVertices[j].x;
  double  yc = fVertices[j].y;
  double  zab = 2 * fDz;
  double  zac = 0;

  if ((std::fabs(xa - xc) + std::fabs(ya - yc)) < halfCarTolerance)
  {
    xc = fVertices[j + 4].x;
    yc = fVertices[j + 4].y;
    zac = 2 * fDz;
    zab = 2 * fDz;

    //Line case
    //
    if ((std::fabs(xb - xc) + std::fabs(yb - yc)) < halfCarTolerance)
    {
      return UUtils::kInfinity;
    }
  }
  double  a = (yb - ya) * zac - (yc - ya) * zab;
  double  b = (xc - xa) * zab - (xb - xa) * zac;
  double  c = (xb - xa) * (yc - ya) - (xc - xa) * (yb - ya);
  double  d = -xa * a - ya * b + fDz * c;
  double  t = a * v.x() + b * v.y() + c * v.z();

  if (t != 0)
  {
    t = -(a * p.x() + b * p.y() + c * p.z() + d) / t;
  }
  if ((t < halfCarTolerance) && (t > -halfCarTolerance))
  {
    if (NormalToPlane(p, ipl).Dot(v) < VUSolid::fgTolerance)
    {
      t = UUtils::kInfinity;
    }
    else
    {
      t = 0;
    }
  }
  if (Inside(p + v * t) != eSurface)
  {
    t = UUtils::kInfinity;
  }

  return t;
}

// --------------------------------------------------------------------
double UGenericTrap::DistanceToOut(const UVector3& p, const UVector3&  v,
                                   UVector3& aNormalVector, bool& aConvex, double) const

{
#ifdef G4TESS_TEST
  if (fTessellatedSolid)
  {
    //return fTessellatedSolid->DistanceToOut(p, v, calcNorm, validNorm, n);
  }
#endif
  static const double  halfCarTolerance = VUSolid::fgTolerance * 0.5;
  double  distmin;
  bool lateral_cross = false;
  ESide side = kUndefined;
  aConvex = true;
  //   *validNorm=true;  // All normals are valid

  if (v.z() < 0)
  {
    distmin = (-fDz - p.z()) / v.z();
    side = kMZ;
    aNormalVector.Set(0, 0, -1);
  }
  else
  {
    if (v.z() > 0)
    {
      distmin = (fDz - p.z()) / v.z();
      side = kPZ;
      aNormalVector.Set(0, 0, 1);
    }
    else
    {
      distmin = UUtils::kInfinity;
    }
  }

  double  dz2 = 0.5 / fDz;
  double  xa, xb, xc, xd;
  double  ya, yb, yc, yd;

  for (int  ipl = 0; ipl < 4; ipl++)
  {
    int  j = (ipl + 1) % 4;
    xa = fVertices[ipl].x;
    ya = fVertices[ipl].y;
    xb = fVertices[ipl + 4].x;
    yb = fVertices[ipl + 4].y;
    xc = fVertices[j].x;
    yc = fVertices[j].y;
    xd = fVertices[4 + j].x;
    yd = fVertices[4 + j].y;

    if (((std::fabs(xb - xd) + std::fabs(yb - yd)) < halfCarTolerance)
        || ((std::fabs(xa - xc) + std::fabs(ya - yc)) < halfCarTolerance))
    {
      double  q = DistToTriangle(p, v, ipl) ;
      if ((q >= 0) && (q < distmin))
      {
        distmin = q;
        lateral_cross = true;
        side = ESide(ipl + 1);
      }
      continue;
    }
    double  tx1 = dz2 * (xb - xa);
    double  ty1 = dz2 * (yb - ya);
    double  tx2 = dz2 * (xd - xc);
    double  ty2 = dz2 * (yd - yc);
    double  dzp = fDz + p.z();
    double  xs1 = xa + tx1 * dzp;
    double  ys1 = ya + ty1 * dzp;
    double  xs2 = xc + tx2 * dzp;
    double  ys2 = yc + ty2 * dzp;
    double  dxs = xs2 - xs1;
    double  dys = ys2 - ys1;
    double  dtx = tx2 - tx1;
    double  dty = ty2 - ty1;
    double  a = (dtx * v.y() - dty * v.x() + (tx1 * ty2 - tx2 * ty1) * v.z()) * v.z();
    double  b = dxs * v.y() - dys * v.x() + (dtx * p.y() - dty * p.x() + ty2 * xs1 - ty1 * xs2
                                         + tx1 * ys2 - tx2 * ys1) * v.z();
    double  c = dxs * p.y() - dys * p.x() + xs1 * ys2 - xs2 * ys1;
    double  q = UUtils::kInfinity;

    if (std::fabs(a) < VUSolid::fgTolerance)
    {
      if (std::fabs(b) < VUSolid::fgTolerance)
      {
        continue;
      }
      q = -c / b;

      // Check for Point on the Surface
      //
      if ((q > -halfCarTolerance) && (q < distmin))
      {
        if (q < halfCarTolerance)
        {
          if (NormalToPlane(p, ipl).Dot(v) < 0.)
          {
            continue;
          }
        }
        distmin = q;
        lateral_cross = true;
        side = ESide(ipl + 1);
      }
      continue;
    }
    double  d = b * b - 4 * a * c;
    if (d >= 0.)
    {
      if (a > 0)
      {
        q = 0.5 * (-b - std::sqrt(d)) / a;
      }
      else
      {
        q = 0.5 * (-b + std::sqrt(d)) / a;
      }

      // Check for Point on the Surface
      //
      if (q > -halfCarTolerance)
      {
        if (q < distmin)
        {
          if (q < halfCarTolerance)
          {
            if (NormalToPlane(p, ipl).Dot(v) < 0.) // Check second root
            {
              if (a > 0)
              {
                q = 0.5 * (-b + std::sqrt(d)) / a;
              }
              else
              {
                q = 0.5 * (-b - std::sqrt(d)) / a;
              }
              if ((q > halfCarTolerance) && (q < distmin))
              {
                distmin = q;
                lateral_cross = true;
                side = ESide(ipl + 1);
              }
              continue;
            }
          }
          distmin = q;
          lateral_cross = true;
          side = ESide(ipl + 1);
        }
      }
      else
      {
        if (a > 0)
        {
          q = 0.5 * (-b + std::sqrt(d)) / a;
        }
        else
        {
          q = 0.5 * (-b - std::sqrt(d)) / a;
        }

        // Check for Point on the Surface
        //
        if ((q > -halfCarTolerance) && (q < distmin))
        {
          if (q < halfCarTolerance)
          {
            if (NormalToPlane(p, ipl).Dot(v) < 0.) // Check second root
            {
              if (a > 0)
              {
                q = 0.5 * (-b - std::sqrt(d)) / a;
              }
              else
              {
                q = 0.5 * (-b + std::sqrt(d)) / a;
              }
              if ((q > halfCarTolerance) && (q < distmin))
              {
                distmin = q;
                lateral_cross = true;
                side = ESide(ipl + 1);
              }
              continue;
            }
          }
          distmin = q;
          lateral_cross = true;
          side = ESide(ipl + 1);
        }
      }
    }
  }
  if (!lateral_cross)  // Make sure that track crosses the top or bottom
  {
    if (distmin >= UUtils::kInfinity)
    {
      distmin = VUSolid::fgTolerance;
    }
    UVector3 pt = p + distmin * v;

    // Check if propagated point is in the polygon
    //
    int  i = 0;
    if (v.z() > 0.)
    {
      i = 4;
    }
    //std::vector<UVector2> xy;
    //for ( int  j=0; j<4; j++)  { xy.push_back(fVertices[i+j]); }
    UVector2 xy[4];
    for (int  j = 0; j < 4; j++)
    {
      xy[j] = fVertices[i + j];
    }

    // Check Inside
    //
    if (InsidePolygone(pt, xy) == eOutside)
    {

      if (v.z() > 0)
      {
        side = kPZ;
        aNormalVector.Set(0, 0, 1);
      }
      else
      {
        side = kMZ;
        aNormalVector.Set(0, 0, -1);
      }

      return 0.;
    }
    else
    {
      if (v.z() > 0)
      {
        side = kPZ;
      }
      else
      {
        side = kMZ;
      }
    }
  }


  UVector3 pt = p + v * distmin;
  switch (side)
  {
    case kXY0:
      aNormalVector = NormalToPlane(pt, 0);
      break;
    case kXY1:
      aNormalVector = NormalToPlane(pt, 1);
      break;
    case kXY2:
      aNormalVector = NormalToPlane(pt, 2);
      break;
    case kXY3:
      aNormalVector = NormalToPlane(pt, 3);
      break;
    case kPZ:
      aNormalVector.Set(0, 0, 1);
      break;
    case kMZ:
      aNormalVector.Set(0, 0, -1);
      break;
    default:
      // DumpInfo();
      std::ostringstream message;
      int  oldprc = message.precision(16);
      message << "Undefined side for valid surface normal to solid." << std::endl
              << "Position:" << std::endl
              << "  p.x() = "   << p.x() << " mm" << std::endl
              << "  p.y() = "   << p.y() << " mm" << std::endl
              << "  p.z() = "   << p.z() << " mm" << std::endl
              << "Direction:" << std::endl
              << "  v.x() = "   << v.x() << std::endl
              << "  v.y() = "   << v.y() << std::endl
              << "  v.z() = "   << v.z() << std::endl
              << "Proposed distance :" << std::endl
              << "  distmin = " << distmin << " mm";
      message.precision(oldprc);
      UUtils::Exception("UGenericTrap::DistanceToOut(p,v,..)",
                        "GeomSolids1002", UWarning, 1, message.str().c_str());
      break;

  }
  if (distmin < halfCarTolerance)
  {
    distmin = 0.;
  }
  //std::cout<<"aNormalVector="<<aNormalVector<<std::endl;
  return distmin;
}
// --------------------------------------------------------------------

double  UGenericTrap::SafetyFromInside(const UVector3& p, bool) const
{

#ifdef G4TESS_TEST
  if (fTessellatedSolid)
  {
    return fTessellatedSolid->DistanceToOut(p);
  }
#endif

  double  safz = fDz - std::fabs(p.z());
  if (safz < 0)
  {
    safz = 0;
  }

  double  safe  = safz;
  double  safxy = safz;

  for (int  iseg = 0; iseg < 4; iseg++)
  {
    safxy = std::fabs(SafetyToFace(p, iseg));
    if (safxy < safe)
    {
      safe = safxy;
    }
  }

  return safe;
}

// --------------------------------------------------------------------
void UGenericTrap::Extent(UVector3& aMin, UVector3& aMax) const
{
#ifdef G4TESS_TEST
  if (fTessellatedSolid)
  {
    //return fTessellatedSolid->CalculateExtent(pAxis, pVoxelLimit,
    pTransform, pMin, pMax);
  }
#endif

  // Computes bounding vectors for a shape
  //

  aMin = GetMinimumBBox();
  aMax = GetMaximumBBox();

}

// --------------------------------------------------------------------

VUSolid* UGenericTrap::Clone() const
{
  return new UGenericTrap(*this);
}

// --------------------------------------------------------------------

std::ostream& UGenericTrap::StreamInfo(std::ostream& os) const
{
  int  oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " *** \n"
     << "    =================================================== \n"
     << " Solid geometry type: " << GetEntityType()  << std::endl
     << "   half length Z: " << fDz << " mm \n"
     << "   list of vertices:\n";

  for (int  i = 0; i < fgkNofVertices; ++i)
  {
    os << std::setw(5) << "#" << i
       << "   vx = " << fVertices[i].x << " mm"
       << "   vy = " << fVertices[i].y << " mm" << std::endl;
  }
  os.precision(oldprc);

  return os;
}

// --------------------------------------------------------------------

UVector3 UGenericTrap::GetPointOnSurface() const
{

#ifdef G4TESS_TEST
  if (fTessellatedSolid)
  {
    //return fTessellatedSolid->GetPointOnSurface();
  }
#endif

  UVector3 point;
  UVector2 u, v, w;
  double  rand, area, chose, cf, lambda0, lambda1, alfa, beta, zp;
  int  ipl, j;

  std::vector<UVector3> vertices;
  for (int  i = 0; i < 4; i++)
  {
    vertices.push_back(UVector3(fVertices[i].x, fVertices[i].y, -fDz));
  }
  for (int  i = 4; i < 8; i++)
  {
    vertices.push_back(UVector3(fVertices[i].x, fVertices[i].y, fDz));
  }

  // Surface Area of Planes(only estimation for twisted)
  //
  double  Surface0 = GetFaceSurfaceArea(vertices[0], vertices[1],
                                        vertices[2], vertices[3]); //-fDz plane
  double  Surface1 = GetFaceSurfaceArea(vertices[0], vertices[1],
                                        vertices[5], vertices[4]); // Lat plane
  double  Surface2 = GetFaceSurfaceArea(vertices[3], vertices[0],
                                        vertices[4], vertices[7]); // Lat plane
  double  Surface3 = GetFaceSurfaceArea(vertices[2], vertices[3],
                                        vertices[7], vertices[6]); // Lat plane
  double  Surface4 = GetFaceSurfaceArea(vertices[2], vertices[1],
                                        vertices[5], vertices[6]); // Lat plane
  double  Surface5 = GetFaceSurfaceArea(vertices[4], vertices[5],
                                        vertices[6], vertices[7]); // fDz plane
  rand = UUtils::Random();
  area = Surface0 + Surface1 + Surface2 + Surface3 + Surface4 + Surface5;
  chose = rand * area;

  if ((chose < Surface0)
      || (chose > (Surface0 + Surface1 + Surface2 + Surface3 + Surface4)))
  {
    // fDz or -fDz Plane
    ipl = int (UUtils::Random() * 4);
    j = (ipl + 1) % 4;
    if (chose < Surface0)
    {
      zp = -fDz;
      u = fVertices[ipl];
      v = fVertices[j];
      w = fVertices[(ipl + 3) % 4];
    }
    else
    {
      zp = fDz;
      u = fVertices[ipl + 4];
      v = fVertices[j + 4];
      w = fVertices[(ipl + 3) % 4 + 4];
    }
    alfa = UUtils::Random();
    beta = UUtils::Random();
    lambda1 = alfa * beta;
    lambda0 = alfa - lambda1;
    v = v - u;
    w = w - u;
    v = u + lambda0 * v + lambda1 * w;
  }
  else                                     // Lateral Plane Twisted or Not
  {
    if (chose < Surface0 + Surface1)
    {
      ipl = 0;
    }
    else if (chose < Surface0 + Surface1 + Surface2)
    {
      ipl = 1;
    }
    else if (chose < Surface0 + Surface1 + Surface2 + Surface3)
    {
      ipl = 2;
    }
    else
    {
      ipl = 3;
    }
    j = (ipl + 1) % 4;
    zp = -fDz + UUtils::Random() * 2 * fDz;
    cf = 0.5 * (fDz - zp) / fDz;
    u = fVertices[ipl + 4] + cf * (fVertices[ipl] - fVertices[ipl + 4]);
    v = fVertices[j + 4] + cf * (fVertices[j] - fVertices[j + 4]);
    v = u + (v - u) * UUtils::Random();
  }
  point = UVector3(v.x, v.y, zp);
  //if(Inside(point)!=eSurface){std::cout<<"GenericTrap::GetPointOnSurface-Point is not"<<point<<" ipl="<<ipl<<std::endl;}
  return point;
}

// --------------------------------------------------------------------

double  UGenericTrap::Capacity()
{
  if (fCubicVolume != 0.)
  {
    ;
  }
  else
  {
    fCubicVolume = VUSolid::Capacity();
  }
  return fCubicVolume;
}

// --------------------------------------------------------------------

double  UGenericTrap::SurfaceArea()
{
  if (fSurfaceArea != 0.)
  {
    ;
  }
  else
  {
    std::vector<UVector3> vertices;
    for (int  i = 0; i < 4; i++)
    {
      vertices.push_back(UVector3(fVertices[i].x, fVertices[i].y, -fDz));
    }
    for (int  i = 4; i < 8; i++)
    {
      vertices.push_back(UVector3(fVertices[i].x, fVertices[i].y, fDz));
    }

    // Surface Area of Planes(only estimation for twisted)
    //
    double  fSurface0 = GetFaceSurfaceArea(vertices[0], vertices[1],
                                           vertices[2], vertices[3]); //-fDz plane
    double  fSurface1 = GetFaceSurfaceArea(vertices[0], vertices[1],
                                           vertices[5], vertices[4]); // Lat plane
    double  fSurface2 = GetFaceSurfaceArea(vertices[3], vertices[0],
                                           vertices[4], vertices[7]); // Lat plane
    double  fSurface3 = GetFaceSurfaceArea(vertices[2], vertices[3],
                                           vertices[7], vertices[6]); // Lat plane
    double  fSurface4 = GetFaceSurfaceArea(vertices[2], vertices[1],
                                           vertices[5], vertices[6]); // Lat plane
    double  fSurface5 = GetFaceSurfaceArea(vertices[4], vertices[5],
                                           vertices[6], vertices[7]); // fDz plane

    // Total Surface Area
    //
    if (!fIsTwisted)
    {
      fSurfaceArea = fSurface0 + fSurface1 + fSurface2
                     + fSurface3 + fSurface4 + fSurface5;
    }
    else
    {
      fSurfaceArea = VUSolid::SurfaceArea();
    }
  }
  return fSurfaceArea;
}

// --------------------------------------------------------------------

double  UGenericTrap::GetFaceSurfaceArea(const UVector3& p0,
                                         const UVector3& p1,
                                         const UVector3& p2,
                                         const UVector3& p3) const
{
  // Auxiliary method for Get Surface Area of Face

  double  aOne, aTwo;
  UVector3 t, u, v, w, Area, normal;

  t = p2 - p1;
  u = p0 - p1;
  v = p2 - p3;
  w = p0 - p3;

  Area = w.Cross(v);
  aOne = 0.5 * Area.Mag();

  Area = t.Cross(u);
  aTwo = 0.5 * Area.Mag();

  return aOne + aTwo;
}

// --------------------------------------------------------------------

bool UGenericTrap::ComputeIsTwisted()
{
  // Computes tangents of twist angles (angles between projections on XY plane
  // of corresponding -dz +dz edges).

  bool twisted = false;
  double  dx1, dy1, dx2, dy2;
  int  nv = fgkNofVertices / 2;

  for (int  i = 0; i < 4; i++)
  {
    dx1 = fVertices[(i + 1) % nv].x - fVertices[i].x;
    dy1 = fVertices[(i + 1) % nv].y - fVertices[i].y;
    if ((dx1 == 0) && (dy1 == 0))
    {
      continue;
    }

    dx2 = fVertices[nv + (i + 1) % nv].x - fVertices[nv + i].x;
    dy2 = fVertices[nv + (i + 1) % nv].y - fVertices[nv + i].y;

    if (dx2 == 0 && dy2 == 0)
    {
      continue;
    }
    double  twist_angle = std::fabs(dy1 * dx2 - dx1 * dy2);
    if (twist_angle < fgkTolerance)
    {
      continue;
    }
    twisted = true;
    SetTwistAngle(i, twist_angle);

    // Check on big angles, potentially navigation problem

    twist_angle = std::acos((dx1 * dx2 + dy1 * dy2)
                            / (std::sqrt(dx1 * dx1 + dy1 * dy1)
                               * std::sqrt(dx2 * dx2 + dy2 * dy2)));

    if (std::fabs(twist_angle) > 0.5 * UUtils::kPi + VUSolid::fgTolerance)
    {
      std::ostringstream message;
      message << "Twisted Angle is bigger than 90 degrees - " << GetName()
              << std::endl
              << "     Potential problem of malformed Solid !" << std::endl
              << "     TwistANGLE = " << twist_angle
              << "*rad  for lateral plane N= " << i;
      UUtils::Exception("UGenericTrap::ComputeIsTwisted()", "GeomSolids1002",
                        UWarning, 4, message.str().c_str());
    }
  }

  return twisted;
}

// --------------------------------------------------------------------

bool UGenericTrap::CheckOrder(const std::vector<UVector2>& vertices) const
{
  // Test if the vertices are in a clockwise order, if not reorder them.
  // Also test if they're well defined without crossing opposite segments

  bool clockwise_order = true;
  double  sum1 = 0.;
  double  sum2 = 0.;
  int  j;

  for (int  i = 0; i < 4; i++)
  {
    j = (i + 1) % 4;
    sum1 += vertices[i].x * vertices[j].y - vertices[j].x * vertices[i].y;
    sum2 += vertices[i + 4].x * vertices[j + 4].y
            - vertices[j + 4].x * vertices[i + 4].y;
  }
  if (sum1 * sum2 < -fgkTolerance)
  {
    std::ostringstream message;
    message << "Lower/upper faces defined with opposite clockwise - "
            << GetName();
    UUtils::Exception("UGenericTrap::CheckOrder()", "GeomSolids0002",
                      UFatalError, 1, message.str().c_str());
  }

  if ((sum1 > 0.) || (sum2 > 0.))
  {
    std::ostringstream message;
    message << "Vertices must be defined in clockwise XY planes - "
            << GetName();
    UUtils::Exception("UGenericTrap::CheckOrder()", "GeomSolids1001",
                      UWarning, 1, message.str().c_str());
    clockwise_order = false;
  }

  // Check for illegal crossings
  //
  bool illegal_cross = false;
  illegal_cross = IsSegCrossingZ(vertices[0], vertices[4],
                                 vertices[1], vertices[5]);

  if (!illegal_cross)
  {
    illegal_cross = IsSegCrossingZ(vertices[2], vertices[6],
                                   vertices[3], vertices[7]);
  }
  // +/- dZ planes
  if (!illegal_cross)
  {
    illegal_cross = IsSegCrossing(vertices[0], vertices[1],
                                  vertices[2], vertices[3]);
  }
  if (!illegal_cross)
  {
    illegal_cross = IsSegCrossing(vertices[0], vertices[3],
                                  vertices[1], vertices[2]);
  }
  if (!illegal_cross)
  {
    illegal_cross = IsSegCrossing(vertices[4], vertices[5],
                                  vertices[6], vertices[7]);
  }
  if (!illegal_cross)
  {
    illegal_cross = IsSegCrossing(vertices[4], vertices[7],
                                  vertices[5], vertices[6]);
  }

  if (illegal_cross)
  {
    std::ostringstream message;
    message << "Malformed polygone with opposite sides - " << GetName();
    UUtils::Exception("UGenericTrap::CheckOrderAndSetup()", "GeomSolids0002",
                      UFatalError, 1, message.str().c_str());
  }
  return clockwise_order;
}

// --------------------------------------------------------------------

void UGenericTrap::ReorderVertices(std::vector<UVector3>& vertices) const
{
  // Reorder the vector of vertices

  std::vector<UVector3> oldVertices(vertices);

  for (int  i = 0; i < int (oldVertices.size()); ++i)
  {
    vertices[i] = oldVertices[oldVertices.size() - 1 - i];
  }
}

// --------------------------------------------------------------------

bool
UGenericTrap::IsSegCrossing(const UVector2& a, const UVector2& b,
                            const UVector2& c, const UVector2& d) const
{
  // Check if segments [A,B] and [C,D] are crossing

  bool stand1 = false;
  bool stand2 = false;
  double  dx1, dx2, xm = 0., ym = 0., a1 = 0., a2 = 0., b1 = 0., b2 = 0.;
  dx1 = (b - a).x;
  dx2 = (d - c).x;

  if (std::fabs(dx1) < fgkTolerance)
  {
    stand1 = true;
  }
  if (std::fabs(dx2) < fgkTolerance)
  {
    stand2 = true;
  }
  if (!stand1)
  {
    a1 = (b.x * a.y - a.x * b.y) / dx1;
    b1 = (b - a).y / dx1;
  }
  if (!stand2)
  {
    a2 = (d.x * c.y - c.x * d.y) / dx2;
    b2 = (d - c).y / dx2;
  }
  if (stand1 && stand2)
  {
    // Segments parallel and vertical
    //
    if (std::fabs(a.x - c.x) < fgkTolerance)
    {
      // Check if segments are overlapping
      //
      if (((c.y - a.y) * (c.y - b.y) < -fgkTolerance)
          || ((d.y - a.y) * (d.y - b.y) < -fgkTolerance)
          || ((a.y - c.y) * (a.y - d.y) < -fgkTolerance)
          || ((b.y - c.y) * (b.y - d.y) < -fgkTolerance))
      {
        return true;
      }

      return false;
    }
    // Different x values
    //
    return false;
  }

  if (stand1)    // First segment vertical
  {
    xm = a.x;
    ym = a2 + b2 * xm;
  }
  else
  {
    if (stand2)  // Second segment vertical
    {
      xm = c.x;
      ym = a1 + b1 * xm;
    }
    else  // Normal crossing
    {
      if (std::fabs(b1 - b2) < fgkTolerance)
      {
        // Parallel segments, are they aligned
        //
        if (std::fabs(c.y - (a1 + b1 * c.x)) > fgkTolerance)
        {
          return false;
        }

        // Aligned segments, are they overlapping
        //
        if (((c.x - a.x) * (c.x - b.x) < -fgkTolerance)
            || ((d.x - a.x) * (d.x - b.x) < -fgkTolerance)
            || ((a.x - c.x) * (a.x - d.x) < -fgkTolerance)
            || ((b.x - c.x) * (b.x - d.x) < -fgkTolerance))
        {
          return true;
        }

        return false;
      }
      xm = (a1 - a2) / (b2 - b1);
      ym = (a1 * b2 - a2 * b1) / (b2 - b1);
    }
  }

  // Check if crossing point is both between A,B and C,D
  //
  double  check = (xm - a.x) * (xm - b.x) + (ym - a.y) * (ym - b.y);
  if (check > -fgkTolerance)
  {
    return false;
  }
  check = (xm - c.x) * (xm - d.x) + (ym - c.y) * (ym - d.y);
  if (check > -fgkTolerance)
  {
    return false;
  }

  return true;
}

// --------------------------------------------------------------------

bool
UGenericTrap::IsSegCrossingZ(const UVector2& a, const UVector2& b,
                             const UVector2& c, const UVector2& d) const
{
  // Check if segments [A,B] and [C,D] are crossing when
  // A and C are on -dZ and B and D are on +dZ

  // Calculate the Intersection point between two lines in 3D
  //
  UVector3 temp1, temp2;
  UVector3 v1, v2, p1, p2, p3, p4, dv;
  double  q, det;
  p1 = UVector3(a.x, a.y, -fDz);
  p2 = UVector3(c.x, c.y, -fDz);
  p3 = UVector3(b.x, b.y, fDz);
  p4 = UVector3(d.x, d.y, fDz);
  v1 = p3 - p1;
  v2 = p4 - p2;
  dv = p2 - p1;

  // In case of Collapsed Vertices No crossing
  //
  if ((std::fabs(dv.x()) < VUSolid::fgTolerance) &&
      (std::fabs(dv.y()) < VUSolid::fgTolerance))
  {
    return false;
  }

  if ((std::fabs((p4 - p3).x()) < VUSolid::fgTolerance) &&
      (std::fabs((p4 - p3).y()) < VUSolid::fgTolerance))
  {
    return false;
  }

  // First estimate if Intersection is possible( if det is 0)
  //
  det = dv.x() * v1.y() * v2.z() + dv.y() * v1.z() * v2.x()
        - dv.x() * v1.z() * v2.y() - dv.y() * v1.x() * v2.z();

  if (std::fabs(det) < VUSolid::fgTolerance) //Intersection
  {
    temp1 = v1.Cross(v2);
    temp2 = (p2 - p1).Cross(v2);
    if (temp1.Dot(temp2) < 0)
    {
      return false;  // intersection negative
    }
    q = temp1.Mag();

    if (q < VUSolid::fgTolerance)
    {
      return false;  // parallel lines
    }
    q = ((dv).Cross(v2)).Mag() / q;

    if (q < 1. - VUSolid::fgTolerance)
    {
      return true;
    }
  }
  return false;
}

// --------------------------------------------------------------------

VUFacet*
UGenericTrap::MakeDownFacet(const std::vector<UVector3>& fromVertices,
                            int  ind1, int  ind2, int  ind3) const
{
  // Create a triangular facet from the polygon points given by indices
  // forming the down side ( the normal goes in -z)
  // Do not create facet if 2 vertices are the same

  if ((fromVertices[ind1] == fromVertices[ind2]) ||
      (fromVertices[ind2] == fromVertices[ind3]) ||
      (fromVertices[ind1] == fromVertices[ind3]))
  {
    return 0;
  }

  std::vector<UVector3> vertices;
  vertices.push_back(fromVertices[ind1]);
  vertices.push_back(fromVertices[ind2]);
  vertices.push_back(fromVertices[ind3]);

  // first vertex most left
  //
  UVector3 cross = (vertices[1] - vertices[0]).Cross(vertices[2] - vertices[1]);

  if (cross.z() > 0.0)
  {
    // Should not happen, as vertices should have been reordered at this stage

    std::ostringstream message;
    message << "Vertices in wrong order - " << GetName();
    UUtils::Exception("UGenericTrap::MakeDownFacet", "GeomSolids0002",
                      UFatalError, 1, message.str().c_str());
  }

  return new UTriangularFacet(vertices[0], vertices[1], vertices[2], UABSOLUTE);
}

// --------------------------------------------------------------------

VUFacet*
UGenericTrap::MakeUpFacet(const std::vector<UVector3>& fromVertices,
                          int  ind1, int  ind2, int  ind3) const
{
  // Create a triangular facet from the polygon points given by indices
  // forming the upper side ( z>0 )

  // Do not create facet if 2 vertices are the same
  //
  if ((fromVertices[ind1] == fromVertices[ind2]) ||
      (fromVertices[ind2] == fromVertices[ind3]) ||
      (fromVertices[ind1] == fromVertices[ind3]))
  {
    return 0;
  }

  std::vector<UVector3> vertices;
  vertices.push_back(fromVertices[ind1]);
  vertices.push_back(fromVertices[ind2]);
  vertices.push_back(fromVertices[ind3]);

  // First vertex most left
  //
  UVector3 cross = (vertices[1] - vertices[0]).Cross(vertices[2] - vertices[1]);

  if (cross.z() < 0.0)
  {
    // Should not happen, as vertices should have been reordered at this stage

    std::ostringstream message;
    message << "Vertices in wrong order - " << GetName();
    UUtils::Exception("UGenericTrap::MakeUpFacet", "GeomSolids0002",
                      UFatalError, 1, message.str().c_str());
  }

  return new UTriangularFacet(vertices[0], vertices[1], vertices[2], UABSOLUTE);
}

// --------------------------------------------------------------------

VUFacet*
UGenericTrap::MakeSideFacet(const UVector3& downVertex0,
                            const UVector3& downVertex1,
                            const UVector3& upVertex1,
                            const UVector3& upVertex0) const
{
  // Creates a triangular facet from the polygon points given by indices
  // forming the upper side ( z>0 )

  if ((downVertex0 == downVertex1) && (upVertex0 == upVertex1))
  {
    return 0;
  }

  if (downVertex0 == downVertex1)
  {
    return new UTriangularFacet(downVertex0, upVertex1, upVertex0, UABSOLUTE);
  }

  if (upVertex0 == upVertex1)
  {
    return new UTriangularFacet(downVertex0, downVertex1, upVertex0, UABSOLUTE);
  }

  return new UQuadrangularFacet(downVertex0, downVertex1,
                                upVertex1, upVertex0, UABSOLUTE);
}

// --------------------------------------------------------------------

UTessellatedSolid* UGenericTrap::CreateTessellatedSolid() const
{
  // 3D vertices
  //
  int  nv = fgkNofVertices / 2;
  std::vector<UVector3> downVertices;
  for (int  i = 0; i < nv; i++)
  {
    downVertices.push_back(UVector3(fVertices[i].x,
                                    fVertices[i].y, -fDz));
  }

  std::vector<UVector3> upVertices;
  for (int  i = nv; i < 2 * nv; i++)
  {
    upVertices.push_back(UVector3(fVertices[i].x,
                                  fVertices[i].y, fDz));
  }

  // Reorder vertices if they are not ordered anti-clock wise
  //
  UVector3 cross
    = (downVertices[1] - downVertices[0]).Cross(downVertices[2] - downVertices[1]);
  UVector3 cross1
    = (upVertices[1] - upVertices[0]).Cross(upVertices[2] - upVertices[1]);
  if ((cross.z() > 0.0) || (cross1.z() > 0.0))
  {
    ReorderVertices(downVertices);
    ReorderVertices(upVertices);
  }

  UTessellatedSolid* tessellatedSolid = new UTessellatedSolid(GetName());

  VUFacet* facet = 0;
  facet = MakeDownFacet(downVertices, 0, 1, 2);
  if (facet)
  {
    tessellatedSolid->AddFacet(facet);
  }
  facet = MakeDownFacet(downVertices, 0, 2, 3);
  if (facet)
  {
    tessellatedSolid->AddFacet(facet);
  }
  facet = MakeUpFacet(upVertices, 0, 2, 1);
  if (facet)
  {
    tessellatedSolid->AddFacet(facet);
  }
  facet = MakeUpFacet(upVertices, 0, 3, 2);
  if (facet)
  {
    tessellatedSolid->AddFacet(facet);
  }

  // The quadrangular sides
  //
  for (int  i = 0; i < nv; ++i)
  {
    int  j = (i + 1) % nv;
    facet = MakeSideFacet(downVertices[j], downVertices[i],
                          upVertices[i], upVertices[j]);

    if (facet)
    {
      tessellatedSolid->AddFacet(facet);
    }
  }

  tessellatedSolid->SetSolidClosed(true);

  return tessellatedSolid;
}

// --------------------------------------------------------------------

void UGenericTrap::ComputeBBox()
{
  // Computes bounding box for a shape.

  double  minX, maxX, minY, maxY;
  minX = maxX = fVertices[0].x;
  minY = maxY = fVertices[0].y;

  for (int  i = 1; i < fgkNofVertices; i++)
  {
    if (minX > fVertices[i].x)
    {
      minX = fVertices[i].x;
    }
    if (maxX < fVertices[i].x)
    {
      maxX = fVertices[i].x;
    }
    if (minY > fVertices[i].y)
    {
      minY = fVertices[i].y;
    }
    if (maxY < fVertices[i].y)
    {
      maxY = fVertices[i].y;
    }
  }
  fMinBBoxVector = UVector3(minX, minY, -fDz);
  fMaxBBoxVector = UVector3(maxX, maxY, fDz);
  fBoundBox = new UBox("BBoxGenericTrap", std::max(std::fabs(minX), std::fabs(maxX)), std::max(std::fabs(minY), std::fabs(maxY)), fDz);
}

