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
// G4TriangularFacet implementation
//
// 31.10.2004, P R Truscott, QinetiQ Ltd, UK - Created.
// 01.08.2007, P R Truscott, QinetiQ Ltd, UK
//                  Significant modification to correct for errors and enhance
//                  based on patches/observations kindly provided by Rickard
//                  Holmberg.
// 12.10.2012, M Gayer, CERN
//                  New implementation reducing memory requirements by 50%,
//                  and considerable CPU speedup together with the new
//                  implementation of G4TessellatedSolid.
// 23.02.2016, E Tcherniaev, CERN
//                  Improved test to detect degenerate (too small or
//                  too narrow) triangles.
// --------------------------------------------------------------------

#include "G4TriangularFacet.hh"

#include "Randomize.hh"
#include "G4TessellatedGeometryAlgorithms.hh"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//
// Definition of triangular facet using absolute vectors to fVertices.
// From this for first vector is retained to define the facet location and
// two relative vectors (E0 and E1) define the sides and orientation of 
// the outward surface normal.
//
G4TriangularFacet::G4TriangularFacet (const G4ThreeVector& vt0,
                                      const G4ThreeVector& vt1,
                                      const G4ThreeVector& vt2,
                                            G4FacetVertexType vertexType)
  : G4VFacet()
{
  fVertices = new vector<G4ThreeVector>(3);

  SetVertex(0, vt0);
  if (vertexType == ABSOLUTE)
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

  G4ThreeVector E1xE2 = fE1.cross(fE2);
  fArea = 0.5 * E1xE2.mag();
  for (G4int i = 0; i < 3; ++i) fIndices[i] = -1;

  fIsDefined = true;
  G4double delta = kCarTolerance; // Set tolerance for checking

  // Check length of edges
  //
  G4double leng1 = fE1.mag();
  G4double leng2 = (fE2-fE1).mag();
  G4double leng3 = fE2.mag();
  if (leng1 <= delta || leng2 <= delta || leng3 <= delta) 
  {
    fIsDefined = false;
  }

  // Check min height of triangle
  //
  if (fIsDefined)
  {
    if (2.*fArea/std::max(std::max(leng1,leng2),leng3) <= delta)
    {
      fIsDefined = false;
    } 
  }

  // Define facet
  //
  if (!fIsDefined)
  {
    ostringstream message;
    message << "Facet is too small or too narrow." << G4endl
            << "Triangle area = " << fArea << G4endl
            << "P0 = " << GetVertex(0) << G4endl
            << "P1 = " << GetVertex(1) << G4endl
            << "P2 = " << GetVertex(2) << G4endl
            << "Side1 length (P0->P1) = " << leng1 << G4endl
            << "Side2 length (P1->P2) = " << leng2 << G4endl
            << "Side3 length (P2->P0) = " << leng3;
    G4Exception("G4TriangularFacet::G4TriangularFacet()",
		"GeomSolids1001", JustWarning, message);
    fSurfaceNormal.set(0,0,0);
    fA = fB = fC = 0.0;
    fDet = 0.0;
    fCircumcentre = vt0 + 0.5*fE1 + 0.5*fE2;
    fArea = fRadius = 0.0;
  }
  else
  { 
    fSurfaceNormal = E1xE2.unit();
    fA   = fE1.mag2();
    fB   = fE1.dot(fE2);
    fC   = fE2.mag2();
    fDet = std::fabs(fA*fC - fB*fB);

    fCircumcentre = 
      vt0 + (E1xE2.cross(fE1)*fC + fE2.cross(E1xE2)*fA) / (2.*E1xE2.mag2());
    fRadius = (fCircumcentre - vt0).mag();
  }
}

///////////////////////////////////////////////////////////////////////////////
//
G4TriangularFacet::G4TriangularFacet ()
  : fSqrDist(0.)
{
  fVertices = new vector<G4ThreeVector>(3);
  G4ThreeVector zero(0,0,0);
  SetVertex(0, zero);
  SetVertex(1, zero);
  SetVertex(2, zero);
  for (G4int i = 0; i < 3; ++i) fIndices[i] = -1;
  fIsDefined = false;
  fSurfaceNormal.set(0,0,0);
  fA = fB = fC = 0;
  fE1 = zero;
  fE2 = zero;
  fDet = 0.0;
  fArea = fRadius = 0.0;
}

///////////////////////////////////////////////////////////////////////////////
//
G4TriangularFacet::~G4TriangularFacet ()
{
  SetVertices(nullptr);
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TriangularFacet::CopyFrom (const G4TriangularFacet& rhs)
{
  char *p = (char *) &rhs;
  copy(p, p + sizeof(*this), (char *)this);

  if (fIndices[0] < 0 && fVertices == nullptr)
  {
    fVertices = new vector<G4ThreeVector>(3);
    for (G4int i = 0; i < 3; ++i) (*fVertices)[i] = (*rhs.fVertices)[i];
  }
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TriangularFacet::MoveFrom (G4TriangularFacet& rhs)
{
  fSurfaceNormal = std::move(rhs.fSurfaceNormal);
  fArea = std::move(rhs.fArea);
  fCircumcentre = std::move(rhs.fCircumcentre);
  fRadius = std::move(rhs.fRadius);
  fIndices = std::move(rhs.fIndices);
  fA = std::move(rhs.fA); fB = std::move(rhs.fB); fC = std::move(rhs.fC);
  fDet = std::move(rhs.fDet);
  fSqrDist = std::move(rhs.fSqrDist);
  fE1 = std::move(rhs.fE1); fE2 = std::move(rhs.fE2);
  fIsDefined = std::move(rhs.fIsDefined);
  fVertices = std::move(rhs.fVertices);
  rhs.fVertices = nullptr;
}

///////////////////////////////////////////////////////////////////////////////
//
G4TriangularFacet::G4TriangularFacet (const G4TriangularFacet& rhs)
  : G4VFacet(rhs)
{
  CopyFrom(rhs);
}

///////////////////////////////////////////////////////////////////////////////
//
G4TriangularFacet::G4TriangularFacet (G4TriangularFacet&& rhs)
  : G4VFacet(rhs)
{
  MoveFrom(rhs);
}

///////////////////////////////////////////////////////////////////////////////
//
G4TriangularFacet&
G4TriangularFacet::operator=(const G4TriangularFacet& rhs)
{
  SetVertices(nullptr);

  if (this != &rhs)
  {
    delete fVertices;
    CopyFrom(rhs);
  }

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
G4TriangularFacet&
G4TriangularFacet::operator=(G4TriangularFacet&& rhs)
{
  SetVertices(nullptr);

  if (this != &rhs)
  {
    delete fVertices;
    MoveFrom(rhs);
  }

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetClone
//
// Simple member function to generate fA duplicate of the triangular facet.
//
G4VFacet* G4TriangularFacet::GetClone ()
{
  G4TriangularFacet* fc =
    new G4TriangularFacet (GetVertex(0), GetVertex(1), GetVertex(2), ABSOLUTE);
  return fc;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetFlippedFacet
//
// Member function to generate an identical facet, but with the normal vector
// pointing at 180 degrees.
//
G4TriangularFacet* G4TriangularFacet::GetFlippedFacet ()
{
  G4TriangularFacet* flipped =
    new G4TriangularFacet (GetVertex(0), GetVertex(1), GetVertex(2), ABSOLUTE);
  return flipped;
}

///////////////////////////////////////////////////////////////////////////////
//
// Distance (G4ThreeVector)
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
G4ThreeVector G4TriangularFacet::Distance (const G4ThreeVector& p)
{
  G4ThreeVector D  = GetVertex(0) - p;
  G4double d = fE1.dot(D);
  G4double e = fE2.dot(D);
  G4double f = D.mag2();
  G4double q = fB*e - fC*d;
  G4double t = fB*d - fA*e;
  fSqrDist = 0.;

  if (q+t <= fDet)
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
          if (-d >= fA) {q = 1.0; fSqrDist = fA + 2.0*d + f;}
          else         {q = -d/fA; fSqrDist = d*q + f;}
        }
        else
        {
          q = 0.0;
          if       (e >= 0.0) {t = 0.0; fSqrDist = f;}
          else if (-e >= fC)   {t = 1.0; fSqrDist = fC + 2.0*e + f;}
          else                {t = -e/fC; fSqrDist = e*t + f;}
        }
      }
      else
      {
        //
        // We are in region 3.
        //
        q = 0.0;
        if      (e >= 0.0) {t = 0.0; fSqrDist = f;}
        else if (-e >= fC)  {t = 1.0; fSqrDist = fC + 2.0*e + f;}
        else               {t = -e/fC; fSqrDist = e*t + f;}
      }
    }
    else if (t < 0.0)
    {
      //
      // We are in region 5.
      //
      t = 0.0;
      if      (d >= 0.0) {q = 0.0; fSqrDist = f;}
      else if (-d >= fA)  {q = 1.0; fSqrDist = fA + 2.0*d + f;}
      else               {q = -d/fA; fSqrDist = d*q + f;}
    }
    else
    {
      //
      // We are in region 0.
      //
      G4double dist = fSurfaceNormal.dot(D);
      fSqrDist = dist*dist;
      return fSurfaceNormal*dist;
    }
  }
  else
  {
    if (q < 0.0)
    {
      //
      // We are in region 2.
      //
      G4double tmp0 = fB + d;
      G4double tmp1 = fC + e;
      if (tmp1 > tmp0)
      {
        G4double numer = tmp1 - tmp0;
        G4double denom = fA - 2.0*fB + fC;
        if (numer >= denom) {q = 1.0; t = 0.0; fSqrDist = fA + 2.0*d + f;}
        else
        {
          q       = numer/denom;
          t       = 1.0 - q;
          fSqrDist = q*(fA*q + fB*t +2.0*d) + t*(fB*q + fC*t + 2.0*e) + f;
        }
      }
      else
      {
        q = 0.0;
        if      (tmp1 <= 0.0) {t = 1.0; fSqrDist = fC + 2.0*e + f;}
        else if (e >= 0.0)    {t = 0.0; fSqrDist = f;}
        else                  {t = -e/fC; fSqrDist = e*t + f;}
      }
    }
    else if (t < 0.0)
    {
      //
      // We are in region 6.
      //
      G4double tmp0 = fB + e;
      G4double tmp1 = fA + d;
      if (tmp1 > tmp0)
      {
        G4double numer = tmp1 - tmp0;
        G4double denom = fA - 2.0*fB + fC;
        if (numer >= denom) {t = 1.0; q = 0.0; fSqrDist = fC + 2.0*e + f;}
        else
        {
          t       = numer/denom;
          q       = 1.0 - t;
          fSqrDist = q*(fA*q + fB*t +2.0*d) + t*(fB*q + fC*t + 2.0*e) + f;
        }
      }
      else
      {
        t = 0.0;
        if      (tmp1 <= 0.0) {q = 1.0; fSqrDist = fA + 2.0*d + f;}
        else if (d >= 0.0)    {q = 0.0; fSqrDist = f;}
        else                  {q = -d/fA; fSqrDist = d*q + f;}
      }
    }
    else
      //
      // We are in region 1.
      //
    {
      G4double numer = fC + e - fB - d;
      if (numer <= 0.0)
      {
        q       = 0.0;
        t       = 1.0;
        fSqrDist = fC + 2.0*e + f;
      }
      else
      {
        G4double denom = fA - 2.0*fB + fC;
        if (numer >= denom) {q = 1.0; t = 0.0; fSqrDist = fA + 2.0*d + f;}
        else
        {
          q       = numer/denom;
          t       = 1.0 - q;
          fSqrDist = q*(fA*q + fB*t + 2.0*d) + t*(fB*q + fC*t + 2.0*e) + f;
        }
      }
    }
  } 
  //
  //
  // Do fA check for rounding errors in the distance-squared.  It appears that
  // the conventional methods for calculating fSqrDist breaks down when very
  // near to or at the surface (as required by transport).
  // We'll therefore also use the magnitude-squared of the vector displacement.
  // (Note that I've also tried to get around this problem by using the
  // existing equations for
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
  G4ThreeVector u = D + q*fE1 + t*fE2;
  G4double u2 = u.mag2();
  //
  // The following (part of the roundoff error check) is from Oliver Merle'q
  // updates.
  //
  if (fSqrDist > u2) fSqrDist = u2;

  return u;
}

///////////////////////////////////////////////////////////////////////////////
//
// Distance (G4ThreeVector, G4double)
//
// Determines the closest distance between point p and the facet.  This makes
// use of G4ThreeVector G4TriangularFacet::Distance, which stores the
// square of the distance in variable fSqrDist.  If approximate methods show 
// the distance is to be greater than minDist, then forget about further
// computation and return fA very large number.
//
G4double G4TriangularFacet::Distance (const G4ThreeVector& p,
                                            G4double minDist)
{
  //
  // Start with quicky test to determine if the surface of the sphere enclosing
  // the triangle is any closer to p than minDist.  If not, then don't bother
  // about more accurate test.
  //
  G4double dist = kInfinity;
  if ((p-fCircumcentre).mag()-fRadius < minDist)
  {
    //
    // It's possible that the triangle is closer than minDist,
    // so do more accurate assessment.
    //
    dist = Distance(p).mag();
  }
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
// Distance (G4ThreeVector, G4double, G4bool)
//
// Determine the distance to point p.  kInfinity is returned if either:
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
G4double G4TriangularFacet::Distance (const G4ThreeVector& p,
                                            G4double minDist,
                                      const G4bool outgoing)
{
  //
  // Start with quicky test to determine if the surface of the sphere enclosing
  // the triangle is any closer to p than minDist.  If not, then don't bother
  // about more accurate test.
  //
  G4double dist = kInfinity;
  if ((p-fCircumcentre).mag()-fRadius < minDist)
  {
    //
    // It's possible that the triangle is closer than minDist,
    // so do more accurate assessment.
    //
    G4ThreeVector v  = Distance(p);
    G4double dist1 = sqrt(fSqrDist);
    G4double dir = v.dot(fSurfaceNormal);
    G4bool wrongSide = (dir > 0.0 && !outgoing) || (dir < 0.0 && outgoing);
    if (dist1 <= kCarTolerance)
    {
      //
      // Point p is very close to triangle.  Check if it's on the wrong side,
      // in which case return distance of 0.0 otherwise .
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
G4double G4TriangularFacet::Extent (const G4ThreeVector axis)
{
  G4double ss = GetVertex(0).dot(axis);
  G4double sp = GetVertex(1).dot(axis);
  if (sp > ss) ss = sp;
  sp = GetVertex(2).dot(axis);
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
// G4TessellatedSolid::Distance(p,v) calculation.
//
// This member function is thanks the valuable work of Rickard Holmberg.  PT.
// However, "gotos" are the Work of the Devil have been exorcised with
// extreme prejudice!!
//
// IMPORTANT NOTE:  These calculations are predicated on v being fA unit
// vector.  If G4TessellatedSolid or other classes call this member function
// with |v| != 1 then there will be errors.
//
G4bool G4TriangularFacet::Intersect (const G4ThreeVector& p,
                                     const G4ThreeVector& v,
                                           G4bool outgoing,
                                           G4double& distance,
                                           G4double& distFromSurface,
                                           G4ThreeVector& normal)
{
  //
  // Check whether the direction of the facet is consistent with the vector v
  // and the need to be outgoing or ingoing.  If inconsistent, disregard and
  // return false.
  //
  G4double w = v.dot(fSurfaceNormal);
  if ((outgoing && w < -dirTolerance) || (!outgoing && w > dirTolerance))
  {
    distance = kInfinity;
    distFromSurface = kInfinity;
    normal.set(0,0,0);
    return false;
  } 
  //
  // Calculate the orthogonal distance from p to the surface containing the
  // triangle.  Then determine if we're on the right or wrong side of the
  // surface (at fA distance greater than kCarTolerance to be consistent with
  // "outgoing".
  //
  const G4ThreeVector& p0 = GetVertex(0);
  G4ThreeVector D  = p0 - p;
  distFromSurface  = D.dot(fSurfaceNormal);
  G4bool wrongSide = (outgoing && distFromSurface < -0.5*kCarTolerance) ||
    (!outgoing && distFromSurface >  0.5*kCarTolerance);
 
  if (wrongSide)
  {
    distance = kInfinity;
    distFromSurface = kInfinity;
    normal.set(0,0,0);
    return false;
  }

  wrongSide = (outgoing && distFromSurface < 0.0)
           || (!outgoing && distFromSurface > 0.0);
  if (wrongSide)
  {
    //
    // We're slightly on the wrong side of the surface.  Check if we're close
    // enough using fA precise distance calculation.
    //
    G4ThreeVector u = Distance(p);
    if (fSqrDist <= kCarTolerance*kCarTolerance)
    {
      //
      // We're very close.  Therefore return fA small negative number
      // to pretend we intersect.
      //
      // distance = -0.5*kCarTolerance
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
      distance = kInfinity;
      distFromSurface = kInfinity;
      normal.set(0,0,0);
      return false;
    }
  }
  if (w < dirTolerance && w > -dirTolerance)
  {
    //
    // The ray is within the plane of the triangle. Project the problem into 2D
    // in the plane of the triangle. First try to create orthogonal unit vectors
    // mu and nu, where mu is fE1/|fE1|.  This is kinda like
    // the original algorithm due to Rickard Holmberg, but with better
    // mathematical justification than the original method ... however,
    // beware Rickard's was less time-consuming.
    //
    // Note that vprime is not fA unit vector.  We need to keep it unnormalised
    // since the values of distance along vprime (s0 and s1) for intersection
    // with the triangle will be used to determine if we cut the plane at the
    // same time.
    //
    G4ThreeVector mu = fE1.unit();
    G4ThreeVector nu = fSurfaceNormal.cross(mu);
    G4TwoVector pprime(p.dot(mu), p.dot(nu));
    G4TwoVector vprime(v.dot(mu), v.dot(nu));
    G4TwoVector P0prime(p0.dot(mu), p0.dot(nu));
    G4TwoVector E0prime(fE1.mag(), 0.0);
    G4TwoVector E1prime(fE2.dot(mu), fE2.dot(nu));
    G4TwoVector loc[2];
    if (G4TessellatedGeometryAlgorithms::IntersectLineAndTriangle2D(pprime,
                                    vprime, P0prime, E0prime, E1prime, loc))
    {
      //
      // There is an intersection between the line and triangle in 2D.
      // Now check which part of the line intersects with the plane
      // containing the triangle in 3D.
      //
      G4double vprimemag = vprime.mag();
      G4double s0        = (loc[0] - pprime).mag()/vprimemag;
      G4double s1        = (loc[1] - pprime).mag()/vprimemag;
      G4double normDist0 = fSurfaceNormal.dot(s0*v) - distFromSurface;
      G4double normDist1 = fSurfaceNormal.dot(s1*v) - distFromSurface;

      if ((normDist0 < 0.0 && normDist1 < 0.0)
       || (normDist0 > 0.0 && normDist1 > 0.0)
       || (normDist0 == 0.0 && normDist1 == 0.0) ) 
      {
        distance        = kInfinity;
        distFromSurface = kInfinity;
        normal.set(0,0,0);
        return false;
      }
      else
      {
        G4double dnormDist = normDist1 - normDist0;
        if (fabs(dnormDist) < DBL_EPSILON)
        {
          distance = s0;
          normal   = fSurfaceNormal;
          if (!outgoing) distFromSurface = -distFromSurface;
          return true;
        }
        else
        {
          distance = s0 - normDist0*(s1-s0)/dnormDist;
          normal   = fSurfaceNormal;
          if (!outgoing) distFromSurface = -distFromSurface;
          return true;
        }
      }
    }
    else
    {
      distance = kInfinity;
      distFromSurface = kInfinity;
      normal.set(0,0,0);
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
  G4ThreeVector pp = p + v*distance;
  G4ThreeVector DD = p0 - pp;
  G4double d = fE1.dot(DD);
  G4double e = fE2.dot(DD);
  G4double ss = fB*e - fC*d;
  G4double t = fB*d - fA*e;

  G4double sTolerance = (fabs(fB)+ fabs(fC) + fabs(d) + fabs(e))*kCarTolerance;
  G4double tTolerance = (fabs(fA)+ fabs(fB) + fabs(d) + fabs(e))*kCarTolerance;
  G4double detTolerance = (fabs(fA)+ fabs(fC) + 2*fabs(fB) )*kCarTolerance;

  //if (ss < 0.0 || t < 0.0 || ss+t > fDet)
  if (ss < -sTolerance || t < -tTolerance || ( ss+t - fDet ) > detTolerance)
  {
    //
    // The intersection is outside of the triangle.
    //
    distance = distFromSurface = kInfinity;
    normal.set(0,0,0);
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
//
// Auxiliary method, returns a uniform random point on the facet
//
G4ThreeVector G4TriangularFacet::GetPointOnFace() const
{
  G4double u = G4UniformRand();
  G4double v = G4UniformRand();
  if (u+v > 1.) { u = 1. - u; v = 1. - v; }
  return GetVertex(0) + u*fE1 + v*fE2;
}

////////////////////////////////////////////////////////////////////////
//
// GetArea
//
// Auxiliary method for returning the surface fArea
//
G4double G4TriangularFacet::GetArea() const
{
  return fArea;
}

////////////////////////////////////////////////////////////////////////
//
G4GeometryType G4TriangularFacet::GetEntityType () const
{
  return "G4TriangularFacet";
}

////////////////////////////////////////////////////////////////////////
//
G4ThreeVector G4TriangularFacet::GetSurfaceNormal () const
{
  return fSurfaceNormal;
}

////////////////////////////////////////////////////////////////////////
//
void G4TriangularFacet::SetSurfaceNormal (G4ThreeVector normal)
{
  fSurfaceNormal = normal;
}
