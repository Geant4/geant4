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
// UTriangularFacet
//
// 22.08.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include <sstream>
#include <algorithm>

#include "UVector2.hh"
#include "UTessellatedSolid.hh"
#include "VUSolid.hh"
#include "UUtils.hh"
#include "UTriangularFacet.hh"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//
// Definition of triangular facet using absolute vectors to vertices.
// From this for first vector is retained to define the facet location and
// two relative vectors (E0 and E1) define the sides and orientation of
// the outward surface normal.
//
UTriangularFacet::UTriangularFacet(const UVector3& vt0, const UVector3& vt1, const UVector3& vt2, UFacetVertexType vertexType)
  : fSqrDist(0.)
{
  fVertices = new vector<UVector3>(3);

  SetVertex(0, vt0);
  if (vertexType == UABSOLUTE)
  {
    SetVertex(1, vt1);
    SetVertex(2, vt2);
    fE1 = vt1 - vt0;
    fE2 = vt2 - vt0;
  }
  else
  {
    SetVertex(1, vt0 + vt1);
    SetVertex(2, vt0 + vt2);
    fE1 = vt1;
    fE2 = vt2;
  }

  for (int i = 0; i < 3; ++i) fIndices[i] = -1;

  double eMag1 = fE1.Mag();
  double eMag2 = fE2.Mag();
  double eMag3 = (fE2 - fE1).Mag();

  if (eMag1 <= kCarTolerance || eMag2 <= kCarTolerance || eMag3 <= kCarTolerance)
  {
    ostringstream message;
    message << "Length of sides of facet are too small." << endl
            << "P[0] = " << GetVertex(0) << endl
            << "P[1] = " << GetVertex(1) << endl
            << "P[2] = " << GetVertex(2) << endl
            << "Side lengths = P[0]->P[1]" << eMag1 << endl
            << "Side lengths = P[0]->P[2]" << eMag2 << endl
            << "Side lengths = P[1]->P[2]" << eMag3;
    UUtils::Exception("UTriangularFacet::UTriangularFacet()", "GeomSolids1001",
                      UWarning, 1, message.str().c_str());
    fIsDefined     = false;
    fSurfaceNormal.Set(0);
    fA = fB = fC = 0.0;
    fDet = 0.0;
    fArea = fRadius = 0.0;
  }
  else
  {
    fIsDefined     = true;
    fSurfaceNormal = fE1.Cross(fE2).Unit();
    fA   = fE1.Mag2();
    fB   = fE1.Dot(fE2);
    fC   = fE2.Mag2();
    fDet = fabs(fA * fC - fB * fB);

    //    sMin = -0.5*kCarTolerance/sqrt(fA);
    //    sMax = 1.0 - sMin;
    //    tMin = -0.5*kCarTolerance/sqrt(fC);
    //    UVector3 vtmp = 0.25 * (fE1 + fE2);

    fArea = 0.5 * (fE1.Cross(fE2)).Mag();
    double lambda0 , lambda1 ;
    if(std::fabs(fArea) < kCarTolerance*kCarTolerance)
    {
      ostringstream message;
      message << "Area of Facet is too small, possible flat triangle!" << endl
              << "  fVertices[0] = " << GetVertex(0) << endl
              << "  fVertices[1] = " << GetVertex(1) << endl
              << "  fVertices[2] = " << GetVertex(2) << endl
              << "Area = " << fArea;
      UUtils::Exception("UTriangularFacet::UTriangularFacet()",
                        "GeomSolids1001", UWarning, 1, message.str().c_str());
      lambda0 = 0.5;
      lambda1 = 0.5;  
    }
    else
    {
      lambda0 = (fA-fB) * fC / (8.0*fArea*fArea);
      lambda1 = (fC-fB) * fA / (8.0*fArea*fArea);
    }
   
    UVector3 p0 = GetVertex(0);
    fCircumcentre = p0 + lambda0 * fE1 + lambda1 * fE2;
    double radiusSqr = (fCircumcentre - p0).Mag2();
    fRadius = sqrt(radiusSqr);
  }
}


UTriangularFacet::UTriangularFacet()
  : fSqrDist(0.)
{
  fVertices = new vector<UVector3>(3);
  UVector3 zero(0, 0, 0);
  SetVertex(0, zero);
  SetVertex(1, zero);
  SetVertex(2, zero);
  for (int i = 0; i < 3; ++i) fIndices[i] = -1;
  fIsDefined = false;
  fSurfaceNormal.Set(0);
  fA = fB = fC = 0;
  fE1 = zero;
  fE2 = zero;
  fDet = 0.0;
  fArea = fRadius = 0;
}

UTriangularFacet::~UTriangularFacet()
{
  SetVertices(NULL);
}


void UTriangularFacet::CopyFrom(const UTriangularFacet& rhs)
{
  char* p = (char*) &rhs;
  copy(p, p + sizeof(*this), (char*)this);

  if (fIndices[0] < 0 && fVertices)
  {
    fVertices = new vector<UVector3>(3);
    for (int i = 0; i < 3; ++i)(*fVertices)[i] = (*rhs.fVertices)[i];
  }
  fIsDefined = rhs.fIsDefined;
  fSurfaceNormal = rhs.fSurfaceNormal;
  fA = rhs.fA; fB = rhs.fB; fC = rhs.fC;
  fE1 = rhs.fE1;
  fE2 = rhs.fE2;
  fDet = rhs.fDet;
  fArea = rhs.fArea; 
  fRadius = rhs.fRadius;
  fSqrDist = rhs.fSqrDist;
}

UTriangularFacet::UTriangularFacet(const UTriangularFacet& rhs) : VUFacet(rhs)
{
  CopyFrom(rhs);
}

UTriangularFacet& UTriangularFacet::operator=(const UTriangularFacet& rhs)
{
  SetVertices(NULL);

  if (this != &rhs)
    CopyFrom(rhs);

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetClone
//
// Simple member function to generate fA duplicate of the triangular facet.
//
VUFacet* UTriangularFacet::GetClone()
{
  UTriangularFacet* fc = new UTriangularFacet(GetVertex(0), GetVertex(1), GetVertex(2), UABSOLUTE);
  return fc;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetFlippedFacet
//
// Member function to generate an identical facet, but with the normal vector
// pointing at 180 degrees.
//
UTriangularFacet* UTriangularFacet::GetFlippedFacet()
{
  UTriangularFacet* flipped = new UTriangularFacet(GetVertex(0), GetVertex(1), GetVertex(2), UABSOLUTE);
  return flipped;
}

///////////////////////////////////////////////////////////////////////////////
//
// Distance (UVector3)
//
// Determines the vector between p and the closest point on the facet to p.
// This is based on the algorithm published in "Geometric Tools for Computer
// Graphics," Philip J Scheider and David H Eberly, Elsevier Science (USA),
// 2003.  at the time of writing, the algorithm is also available in fA
// technical note "Distance between point and triangle in 3D," by David Eberly
// at http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
//
// The by-product is the square-distance fSqrDist, which is retained
// in case needed by the other "Distance" member functions.
//
UVector3 UTriangularFacet::Distance(const UVector3& p)
{
  UVector3 D  = GetVertex(0) - p;
  double d = fE1.Dot(D);
  double e = fE2.Dot(D);
  double f = D.Mag2();
  double q = fB * e - fC * d;
  double t = fB * d - fA * e;
  fSqrDist = 0.;

  if (q + t <= fDet)
  {
    if (q < 0.0)
    {
      if (t < 0.0)
      {
        //
        // We are in region 4.
        //
        if (d < 0.0)
        {
          t = 0.0;
          if (-d >= fA)
          {
            q = 1.0;
            fSqrDist = fA + 2.0 * d + f;
          }
          else
          {
            q = -d / fA;
            fSqrDist = d * q + f;
          }
        }
        else
        {
          q = 0.0;
          if (e >= 0.0)
          {
            t = 0.0;
            fSqrDist = f;
          }
          else if (-e >= fC)
          {
            t = 1.0;
            fSqrDist = fC + 2.0 * e + f;
          }
          else
          {
            t = -e / fC;
            fSqrDist = e * t + f;
          }
        }
      }
      else
      {
        //
        // We are in region 3.
        //
        q = 0.0;
        if (e >= 0.0)
        {
          t = 0.0;
          fSqrDist = f;
        }
        else if (-e >= fC)
        {
          t = 1.0;
          fSqrDist = fC + 2.0 * e + f;
        }
        else
        {
          t = -e / fC;
          fSqrDist = e * t + f;
        }
      }
    }
    else if (t < 0.0)
    {
      //
      // We are in region 5.
      //
      t = 0.0;
      if (d >= 0.0)
      {
        q = 0.0;
        fSqrDist = f;
      }
      else if (-d >= fA)
      {
        q = 1.0;
        fSqrDist = fA + 2.0 * d + f;
      }
      else
      {
        q = -d / fA;
        fSqrDist = d * q + f;
      }
    }
    else
    {
      //
      // We are in region 0.
      //
      q       = q / fDet;
      t       = t / fDet;
      fSqrDist = q * (fA * q + fB * t + 2.0 * d) + t * (fB * q + fC * t + 2.0 * e) + f;
    }
  }
  else
  {
    if (q < 0.0)
    {
      //
      // We are in region 2.
      //
      double tmp0 = fB + d;
      double tmp1 = fC + e;
      if (tmp1 > tmp0)
      {
        double numer = tmp1 - tmp0;
        double denom = fA - 2.0 * fB + fC;
        if (numer >= denom)
        {
          q = 1.0;
          t = 0.0;
          fSqrDist = fA + 2.0 * d + f;
        }
        else
        {
          q       = numer / denom;
          t       = 1.0 - q;
          fSqrDist = q * (fA * q + fB * t + 2.0 * d) + t * (fB * q + fC * t + 2.0 * e) + f;
        }
      }
      else
      {
        q = 0.0;
        if (tmp1 <= 0.0)
        {
          t = 1.0;
          fSqrDist = fC + 2.0 * e + f;
        }
        else if (e >= 0.0)
        {
          t = 0.0;
          fSqrDist = f;
        }
        else
        {
          t = -e / fC;
          fSqrDist = e * t + f;
        }
      }
    }
    else if (t < 0.0)
    {
      //
      // We are in region 6.
      //
      double tmp0 = fB + e;
      double tmp1 = fA + d;
      if (tmp1 > tmp0)
      {
        double numer = tmp1 - tmp0;
        double denom = fA - 2.0 * fB + fC;
        if (numer >= denom)
        {
          t = 1.0;
          q = 0.0;
          fSqrDist = fC + 2.0 * e + f;
        }
        else
        {
          t       = numer / denom;
          q       = 1.0 - t;
          fSqrDist = q * (fA * q + fB * t + 2.0 * d) + t * (fB * q + fC * t + 2.0 * e) + f;
        }
      }
      else
      {
        t = 0.0;
        if (tmp1 <= 0.0)
        {
          q = 1.0;
          fSqrDist = fA + 2.0 * d + f;
        }
        else if (d >= 0.0)
        {
          q = 0.0;
          fSqrDist = f;
        }
        else
        {
          q = -d / fA;
          fSqrDist = d * q + f;
        }
      }
    }
    else
      //
      // We are in region 1.
      //
    {
      double numer = fC + e - fB - d;
      if (numer <= 0.0)
      {
        q       = 0.0;
        t       = 1.0;
        fSqrDist = fC + 2.0 * e + f;
      }
      else
      {
        double denom = fA - 2.0 * fB + fC;
        if (numer >= denom)
        {
          q = 1.0;
          t = 0.0;
          fSqrDist = fA + 2.0 * d + f;
        }
        else
        {
          q       = numer / denom;
          t       = 1.0 - q;
          fSqrDist = q * (fA * q + fB * t + 2.0 * d) + t * (fB * q + fC * t + 2.0 * e) + f;
        }
      }
    }
  }
  //
  //
  // Do fA check for rounding errors in the distance-squared.  It appears that
  // the conventional methods for calculating fSqrDist breaks down when very near
  // to or at the surface (as required by transport).  We'll therefore also use
  // the magnitude-squared of the vector displacement.  (Note that I've also
  // tried to get around this problem by using the existing equations for
  //
  //    fSqrDist = function(fA,fB,fC,d,q,t)
  //
  // and use fA more accurate addition process which minimises errors and
  // breakdown of cummutitivity [where (A+B)+C != A+(B+C)] but this still
  // doesn't work.
  // Calculation from u = D + q*fE1 + t*fE2 is less efficient, but appears
  // more robust.
  //
  if (fSqrDist < 0.0) fSqrDist = 0.;
  UVector3 u = D + q * fE1 + t * fE2;
  double u2 = u.Mag2();
  //
  // The following (part of the roundoff error check) is from Oliver Merle'q
  // updates.
  //
  if (fSqrDist > u2) fSqrDist = u2;

  return u;
}

///////////////////////////////////////////////////////////////////////////////
//
// Distance (UVector3, double)
//
// Determines the closest distance between point p and the facet.  This makes
// use of UVector3 UTriangularFacet::Distance, which stores the
// square of the distance in variable fSqrDist.  If approximate methods show
// the distance is to be greater than minDist, then forget about further
// computation and return fA very large number.
//
double UTriangularFacet::Distance(const UVector3& p, const double minDist)
{
  //
  // Start with quicky test to determine if the surface of the sphere enclosing
  // the triangle is any closer to p than minDist.  If not, then don't bother
  // about more accurate test.
  //
  double dist = UUtils::kInfinity;
  if ((p - fCircumcentre).Mag() - fRadius < minDist)
  {
    //
    // It's possible that the triangle is closer than minDist, so do more accurate
    // assessment.
    //
    dist = Distance(p).Mag();
    //    dist = sqrt(fSqrDist);
  }
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
// Distance (UVector3, double, bool)
//
// Determine the distance to point p.  UUtils::kInfinity is returned if either:
// (1) outgoing is TRUE and the dot product of the normal vector to the facet
//     and the displacement vector from p to the triangle is negative.
// (2) outgoing is FALSE and the dot product of the normal vector to the facet
//     and the displacement vector from p to the triangle is positive.
// If approximate methods show the distance is to be greater than minDist, then
// forget about further computation and return fA very large number.
//
// This method has been heavily modified thanks to the valuable comments and
// corrections of Rickard Holmberg.
//
double UTriangularFacet::Distance(const UVector3& p, const double minDist, const bool outgoing)
{
  //
  // Start with quicky test to determine if the surface of the sphere enclosing
  // the triangle is any closer to p than minDist.  If not, then don't bother
  // about more accurate test.
  //
  double dist = UUtils::kInfinity;
  if ((p - fCircumcentre).Mag() - fRadius < minDist)
  {
    //
    // It's possible that the triangle is closer than minDist, so do more accurate
    // assessment.
    //
    UVector3 v  = Distance(p);
    double dist1 = sqrt(fSqrDist);
    double dir = v.Dot(fSurfaceNormal);
    bool wrongSide = (dir > 0.0 && !outgoing) || (dir < 0.0 && outgoing);
    if (dist1 <= kCarTolerance)
    {
      //
      // Point p is very close to triangle.  Check if it's on the wrong side, in
      // which case return distance of 0.0 otherwise .
      //
      if (wrongSide) dist = 0.0;
      else dist = dist1;
    }
    else if (!wrongSide) dist = dist1;
  }
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
// Extent
//
// Calculates the furthest the triangle extends in fA particular direction
// defined by the vector axis.
//
double UTriangularFacet::Extent(const UVector3 axis)
{
  double ss = GetVertex(0).Dot(axis);
  double sp = GetVertex(1).Dot(axis);
  if (sp > ss) ss = sp;
  sp = GetVertex(2).Dot(axis);
  if (sp > ss) ss = sp;
  return ss;
}

///////////////////////////////////////////////////////////////////////////////
//
// Intersect
//
// Member function to find the next intersection when going from p in the
// direction of v.  If:
// (1) "outgoing" is TRUE, only consider the face if we are going out through
//     the face.
// (2) "outgoing" is FALSE, only consider the face if we are going in through
//     the face.
// Member functions returns TRUE if there is an intersection, FALSE otherwise.
// Sets the distance (distance along w), distFromSurface (orthogonal distance)
// and normal.
//
// Also considers intersections that happen with negative distance for small
// distances of distFromSurface = 0.5*kCarTolerance in the wrong direction.
// This is to detect kSurface without doing fA full Inside(p) in
// UTessellatedSolid::Distance(p,v) calculation.
//
// This member function is thanks the valuable work of Rickard Holmberg.  PT.
// However, "gotos" are the Work of the Devil have been exorcised with
// extreme prejudice!!
//
// IMPORTANT NOTE:  These calculations are predicated on v being fA unit
// vector.  If UTessellatedSolid or other classes call this member function
// with |v| != 1 then there will be errors.
//
bool UTriangularFacet::Intersect(const UVector3& p, const UVector3& v, bool outgoing, double& distance, double& distFromSurface, UVector3& normal)
{
  //
  // Check whether the direction of the facet is consistent with the vector v
  // and the need to be outgoing or ingoing.  If inconsistent, disregard and
  // return false.
  //
  double w = v.Dot(fSurfaceNormal);
  if ((outgoing && w < -dirTolerance) || (!outgoing && w > dirTolerance))
  {
    distance = UUtils::kInfinity;
    distFromSurface = UUtils::kInfinity;
    normal.Set(0);
    return false;
  }
  //
  // Calculate the orthogonal distance from p to the surface containing the
  // triangle.  Then determine if we're on the right or wrong side of the
  // surface (at fA distance greater than kCarTolerance to be consistent with
  // "outgoing".
  //
  const UVector3& p0 = GetVertex(0);
  UVector3 D  = p0 - p;
  distFromSurface  = D.Dot(fSurfaceNormal);
  bool wrongSide = (outgoing && distFromSurface < -0.5 * kCarTolerance) ||
                   (!outgoing && distFromSurface >  0.5 * kCarTolerance);
  if (wrongSide)
  {
    distance = UUtils::kInfinity;
    distFromSurface = UUtils::kInfinity;
    normal.Set(0);
    return false;
  }

  wrongSide = (outgoing && distFromSurface < 0.0) || (!outgoing && distFromSurface > 0.0);
  if (wrongSide)
  {
    //
    // We're slightly on the wrong side of the surface.  Check if we're close
    // enough using fA precise distance calculation.
    //
    //UVector3 u = Distance(p);
    if (fSqrDist <= kCarTolerance * kCarTolerance)
    {
      //
      // We're very close.  Therefore return fA small negative number to pretend
      // we intersect.
      //
      //      distance = -0.5*kCarTolerance
      distance = 0.0;
      normal = fSurfaceNormal;
      return true;
    }
    else
    {
      //
      // We're close to the surface containing the triangle, but sufficiently
      // far from the triangle, and on the wrong side compared to the directions
      // of the surface normal and v.  There is no intersection.
      //
      distance = UUtils::kInfinity;
      distFromSurface = UUtils::kInfinity;
      normal.Set(0);
      return false;
    }
  }
  if (w < dirTolerance && w > -dirTolerance)
  {
    //
    // The ray is within the plane of the triangle.  Project the problem into 2D
    // in the plane of the triangle.  First try to create orthogonal unit vectors
    // mu and nu, where mu is fE1/|fE1|.  This is kinda like
    // the original algorithm due to Rickard Holmberg, but with better mathematical
    // justification than the original method ... however, beware Rickard's was less
    // time-consuming.
    //
    // Note that vprime is not fA unit vector.  We need to keep it unnormalised
    // since the values of distance along vprime (s0 and s1) for intersection with
    // the triangle will be used to determine if we cut the plane at the same
    // time.
    //
    UVector3 mu = fE1.Unit();
    UVector3 nu = fSurfaceNormal.Cross(mu);
    UVector2 pprime(p.Dot(mu), p.Dot(nu));
    UVector2 vprime(v.Dot(mu), v.Dot(nu));
    UVector2 P0prime(p0.Dot(mu), p0.Dot(nu));
    UVector2 E0prime(fE1.Mag(), 0.0);
    UVector2 E1prime(fE2.Dot(mu), fE2.Dot(nu));
    UVector2 loc[2];
    if (UTessellatedGeometryAlgorithms::IntersectLineAndTriangle2D(pprime, vprime, P0prime, E0prime, E1prime, loc))
    {
      //
      // There is an intersection between the line and triangle in 2D.  Now check
      // which part of the line intersects with the plane containing the triangle
      // in 3D.
      //
      double vprimemag = vprime.mag();
      double s0        = (loc[0] - pprime).mag() / vprimemag;
      double s1        = (loc[1] - pprime).mag() / vprimemag;
      double normDist0 = fSurfaceNormal.Dot(s0 * v) - distFromSurface;
      double normDist1 = fSurfaceNormal.Dot(s1 * v) - distFromSurface;

      if ((normDist0 < 0.0 && normDist1 < 0.0) || (normDist0 > 0.0 && normDist1 > 0.0) ||
          (normDist0 == 0.0 && normDist1 == 0.0))
      {
        distance        = UUtils::kInfinity;
        distFromSurface = UUtils::kInfinity;
        normal.Set(0);
        return false;
      }
      else
      {
        double dnormDist = normDist1 - normDist0;
        if (fabs(dnormDist) < DBL_EPSILON)
        {
          distance = s0;
          normal   = fSurfaceNormal;
          if (!outgoing) distFromSurface = -distFromSurface;
          return true;
        }
        else
        {
          distance = s0 - normDist0 * (s1 - s0) / dnormDist;
          normal   = fSurfaceNormal;
          if (!outgoing) distFromSurface = -distFromSurface;
          return true;
        }
      }

      //      UVector3 dloc   = loc1 - loc0;
      //      UVector3 dlocXv = dloc.cross(v);
      //      double dlocXvmag   = dlocXv.mag();
      //      if (dloc.mag() <= 0.5*kCarTolerance || dlocXvmag <= DBL_EPSILON)
      //      {
      //        distance = loc0.mag();
      //        normal = fSurfaceNormal;
      //        if (!outgoing) distFromSurface = -distFromSurface;
      //        return true;
      //      }

      //      UVector3 loc0Xv   = loc0.cross(v);
      //      UVector3 loc1Xv   = loc1.cross(v);
      //      double sameDir       = -loc0Xv.dot(loc1Xv);
      //      if (sameDir < 0.0)
      //      {
      //        distance        = UUtils::kInfinity;
      //        distFromSurface = UUtils::kInfinity;
      //        normal          = UVector3(0.0,0.0,0.0);
      //        return false;
      //      }
      //      else
      //      {
      //        distance = loc0.mag() + loc0Xv.mag() * dloc.mag()/dlocXvmag;
      //        normal   = fSurfaceNormal;
      //        if (!outgoing) distFromSurface = -distFromSurface;
      //        return true;
      //      }
    }
    else
    {
      distance = UUtils::kInfinity;
      distFromSurface = UUtils::kInfinity;
      normal.Set(0);
      return false;
    }
  }
  //
  //
  // Use conventional algorithm to determine the whether there is an
  // intersection.  This involves determining the point of intersection of the
  // line with the plane containing the triangle, and then calculating if the
  // point is within the triangle.
  //
  distance = distFromSurface / w;
  UVector3 pp = p + v * distance;
  UVector3 DD = p0 - pp;
  double d = fE1.Dot(DD);
  double e = fE2.Dot(DD);
  double ss = fB * e - fC * d;
  double t = fB * d - fA * e;

  double sTolerance = (fabs(fB) + fabs(fC) + fabs(d) + fabs(e)) * kCarTolerance;
  double tTolerance = (fabs(fA) + fabs(fB) + fabs(d) + fabs(e)) * kCarTolerance;
  double detTolerance = (fabs(fA) + fabs(fC) + 2 * fabs(fB)) * kCarTolerance;

  //if (ss < 0.0 || t < 0.0 || ss+t > fDet)
  if (ss < -sTolerance || t < -tTolerance || (ss + t - fDet) > detTolerance)
  {
    //
    // The intersection is outside of the triangle.
    //
    distance = distFromSurface = UUtils::kInfinity;
    normal.Set(0);
    return false;
  }
  else
  {
    //
    // There is an intersection.  Now we only need to set the surface normal.
    //
    normal = fSurfaceNormal;
    if (!outgoing) distFromSurface = -distFromSurface;
    return true;
  }
}

////////////////////////////////////////////////////////////////////////
//
// GetPointOnFace
//```
// Auxiliary method for get fA random point on surface

UVector3 UTriangularFacet::GetPointOnFace() const
{
  double alpha = UUtils::Random(0., 1.);
  double beta = UUtils::Random(0., 1.);
  double lambda1 = alpha * beta;
  double lambda0 = alpha - lambda1;

  return GetVertex(0) + lambda0 * fE1 + lambda1 * fE2;
}

////////////////////////////////////////////////////////////////////////
//
// GetArea
//
// Auxiliary method for returning the surface fArea

double UTriangularFacet::GetArea()
{
  return fArea;
}

UGeometryType UTriangularFacet::GetEntityType() const
{
  return "TriangularFacet";
}

UVector3 UTriangularFacet::GetSurfaceNormal() const
{
  return fSurfaceNormal;
}

void UTriangularFacet::SetSurfaceNormal(UVector3 normal)
{
  fSurfaceNormal = normal;
}
