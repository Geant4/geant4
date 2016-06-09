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
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * subject to DEFCON 705 IPR conditions.                            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4TriangularFacet.cc,v 1.16 2010-09-23 10:27:25 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4TriangularFacet.cc
//
// Date:                15/06/2005
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            UK Ministry of Defence : RAO CRP TD Electronic Systems
// Contract:            C/MAT/N03517
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK - Created.
//
// 01 August 2007   P R Truscott, QinetiQ Ltd, UK
//                  Significant modification to correct for errors and enhance
//                  based on patches/observations kindly provided by Rickard
//                  Holmberg
//
// 26 September 2007
//                  P R Truscott, QinetiQ Ltd, UK
//                  Further chamges implemented to the Intersect member
//                  function to correctly treat rays nearly parallel to the
//                  plane of the triangle.
//
// 12 April 2010    P R Truscott, QinetiQ, bug fixes to treat optical
//                  photon transport, in particular internal reflection
//                  at surface.
//
// 22 August 2011   I Hrivnacova, Orsay, fix in Intersect() to take into
//                  account geometrical tolerance and cases of zero distance
//                  from surface.
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "G4TriangularFacet.hh"
#include "G4TwoVector.hh"
#include "globals.hh"
#include "Randomize.hh"

///////////////////////////////////////////////////////////////////////////////
//
// Definition of triangular facet using absolute vectors to vertices.
// From this for first vector is retained to define the facet location and
// two relative vectors (E0 and E1) define the sides and orientation of 
// the outward surface normal.
//
G4TriangularFacet::G4TriangularFacet (const G4ThreeVector Pt0,
             const G4ThreeVector vt1, const G4ThreeVector vt2,
                   G4FacetVertexType vertexType)
  : G4VFacet(), sMin(0.), sMax(1.), tMin(0.), sqrDist(0.)
{
  tGeomAlg  = G4TessellatedGeometryAlgorithms::GetInstance();
  P0        = Pt0;
  nVertices = 3;
  if (vertexType == ABSOLUTE)
  {
    P.push_back(vt1);
    P.push_back(vt2);
  
    E.push_back(vt1 - P0);
    E.push_back(vt2 - P0);
  }
  else
  {
    P.push_back(P0 + vt1);
    P.push_back(P0 + vt2);
  
    E.push_back(vt1);
    E.push_back(vt2);
  }

  G4double Emag1 = E[0].mag();
  G4double Emag2 = E[1].mag();
  G4double Emag3 = (E[1]-E[0]).mag();
  
  if (Emag1 <= kCarTolerance || Emag2 <= kCarTolerance ||
      Emag3 <= kCarTolerance)
  {
    std::ostringstream message;
    message << "Length of sides of facet are too small." << G4endl
            << "P0 = " << P0   << G4endl
            << "P1 = " << P[0] << G4endl
            << "P2 = " << P[1] << G4endl
            << "Side lengths = P0->P1" << Emag1 << G4endl
            << "Side lengths = P0->P2" << Emag2 << G4endl
            << "Side lengths = P1->P2" << Emag3;
    G4Exception("G4TriangularFacet::G4TriangularFacet()",
                "InvalidSetup", JustWarning, G4String(message.str()));
    isDefined     = false;
    geometryType  = "G4TriangularFacet";
    surfaceNormal = G4ThreeVector(0.0,0.0,0.0);
    a   = 0.0;
    b   = 0.0;
    c   = 0.0;
    det = 0.0;
  }
  else
  {
    isDefined     = true;
    geometryType  = "G4TriangularFacet";
    surfaceNormal = E[0].cross(E[1]).unit();
    a   = E[0].mag2();
    b   = E[0].dot(E[1]);
    c   = E[1].mag2();
    det = std::fabs(a*c - b*b);
    
    sMin = -0.5*kCarTolerance/std::sqrt(a);
    sMax = 1.0 - sMin;
    tMin = -0.5*kCarTolerance/std::sqrt(c);
    
    area = 0.5 * (E[0].cross(E[1])).mag();

//    G4ThreeVector vtmp = 0.25 * (E[0] + E[1]);
    G4double lambda0 = (a-b) * c / (8.0*area*area);
    G4double lambda1 = (c-b) * a / (8.0*area*area);
    circumcentre     = P0 + lambda0*E[0] + lambda1*E[1];
    radiusSqr        = (circumcentre-P0).mag2();
    radius           = std::sqrt(radiusSqr);
  
    for (size_t i=0; i<3; i++) { I.push_back(0); }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// ~G4TriangularFacet
//
// A pretty boring destructor indeed!
//
G4TriangularFacet::~G4TriangularFacet ()
{
  P.clear();
  E.clear();
  I.clear();
}

///////////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4TriangularFacet::G4TriangularFacet (const G4TriangularFacet &rhs)
  : G4VFacet(rhs), a(rhs.a), b(rhs.b), c(rhs.c), det(rhs.det),
    sMin(rhs.sMin), sMax(rhs.sMax), tMin(rhs.tMin), sqrDist(rhs.sqrDist)
{
  tGeomAlg = G4TessellatedGeometryAlgorithms::GetInstance();
}

///////////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
const G4TriangularFacet &G4TriangularFacet::operator=(G4TriangularFacet &rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VFacet::operator=(rhs);

   // Copy data
   //
   a = rhs.a; b = rhs.b; c = rhs.c; det = rhs.det;
   sMin = rhs.sMin; sMax = rhs.sMax; tMin = rhs.tMin; sqrDist = rhs.sqrDist;
   tGeomAlg = G4TessellatedGeometryAlgorithms::GetInstance();

   return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetClone
//
// Simple member function to generate a diplicate of the triangular facet.
//
G4VFacet *G4TriangularFacet::GetClone ()
{
  G4TriangularFacet *fc = new G4TriangularFacet (P0, P[0], P[1], ABSOLUTE);
  G4VFacet *cc         = 0;
  cc                   = fc;
  return cc;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetFlippedFacet
//
// Member function to generate an identical facet, but with the normal vector
// pointing at 180 degrees.
//
G4TriangularFacet *G4TriangularFacet::GetFlippedFacet ()
{
  G4TriangularFacet *flipped = new G4TriangularFacet (P0, P[1], P[0], ABSOLUTE);
  return flipped;
}

///////////////////////////////////////////////////////////////////////////////
//
// Distance (G4ThreeVector)
//
// Determines the vector between p and the closest point on the facet to p.
// This is based on the algorithm published in "Geometric Tools for Computer
// Graphics," Philip J Scheider and David H Eberly, Elsevier Science (USA),
// 2003.  at the time of writing, the algorithm is also available in a
// technical note "Distance between point and triangle in 3D," by David Eberly
// at http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
//
// The by-product is the square-distance sqrDist, which is retained
// in case needed by the other "Distance" member functions.
//
G4ThreeVector G4TriangularFacet::Distance (const G4ThreeVector &p)
{
  G4ThreeVector D  = P0 - p;
  G4double d       = E[0].dot(D);
  G4double e       = E[1].dot(D);
  G4double f       = D.mag2();
  G4double s       = b*e - c*d;
  G4double t       = b*d - a*e;

  sqrDist          = 0.0;

  if (s+t <= det)
  {
    if (s < 0.0)
    {
      if (t < 0.0)
      {
  //
  // We are in region 4.
  //
        if (d < 0.0)
        {
          t = 0.0;
          if (-d >= a) {s = 1.0; sqrDist = a + 2.0*d + f;}
          else         {s = -d/a; sqrDist = d*s + f;}
        }
        else
        {
          s = 0.0;
          if       (e >= 0.0) {t = 0.0; sqrDist = f;}
          else if (-e >= c)   {t = 1.0; sqrDist = c + 2.0*e + f;}
          else                {t = -e/c; sqrDist = e*t + f;}
        }
      }
      else
      {
   //
   // We are in region 3.
   //
        s = 0.0;
        if      (e >= 0.0) {t = 0.0; sqrDist = f;}
        else if (-e >= c)  {t = 1.0; sqrDist = c + 2.0*e + f;}
        else               {t = -e/c; sqrDist = e*t + f;}
      }
    }
    else if (t < 0.0)
    {
   //
   // We are in region 5.
   //
      t = 0.0;
      if      (d >= 0.0) {s = 0.0; sqrDist = f;}
      else if (-d >= a)  {s = 1.0; sqrDist = a + 2.0*d + f;}
      else               {s = -d/a; sqrDist = d*s + f;}
    }
    else
    {
   //
   // We are in region 0.
   //
      s       = s / det;
      t       = t / det;
      sqrDist = s*(a*s + b*t + 2.0*d) + t*(b*s + c*t + 2.0*e) + f;
    }
  }
  else
  {
    if (s < 0.0)
    {
   //
   // We are in region 2.
   //
      G4double tmp0 = b + d;
      G4double tmp1 = c + e;
      if (tmp1 > tmp0)
      {
        G4double numer = tmp1 - tmp0;
        G4double denom = a - 2.0*b + c;
        if (numer >= denom) {s = 1.0; t = 0.0; sqrDist = a + 2.0*d + f;}
        else
        {
          s       = numer/denom;
          t       = 1.0 - s;
          sqrDist = s*(a*s + b*t +2.0*d) + t*(b*s + c*t + 2.0*e) + f;
        }
      }
      else
      {
        s = 0.0;
        if      (tmp1 <= 0.0) {t = 1.0; sqrDist = c + 2.0*e + f;}
        else if (e >= 0.0)    {t = 0.0; sqrDist = f;}
        else                  {t = -e/c; sqrDist = e*t + f;}
      }
    }
    else if (t < 0.0)
    {
   //
   // We are in region 6.
   //
      G4double tmp0 = b + e;
      G4double tmp1 = a + d;
      if (tmp1 > tmp0)
      {
        G4double numer = tmp1 - tmp0;
        G4double denom = a - 2.0*b + c;
        if (numer >= denom) {t = 1.0; s = 0.0; sqrDist = c + 2.0*e + f;}
        else
        {
          t       = numer/denom;
          s       = 1.0 - t;
          sqrDist = s*(a*s + b*t +2.0*d) + t*(b*s + c*t + 2.0*e) + f;
        }
      }
      else
      {
        t = 0.0;
        if      (tmp1 <= 0.0) {s = 1.0; sqrDist = a + 2.0*d + f;}
        else if (d >= 0.0)    {s = 0.0; sqrDist = f;}
        else                  {s = -d/a; sqrDist = d*s + f;}
      }
    }
    else
   //
   // We are in region 1.
   //
    {
      G4double numer = c + e - b - d;
      if (numer <= 0.0)
      {
        s       = 0.0;
        t       = 1.0;
        sqrDist = c + 2.0*e + f;
      }
      else
      {
        G4double denom = a - 2.0*b + c;
        if (numer >= denom) {s = 1.0; t = 0.0; sqrDist = a + 2.0*d + f;}
        else
        {
          s       = numer/denom;
          t       = 1.0 - s;
          sqrDist = s*(a*s + b*t + 2.0*d) + t*(b*s + c*t + 2.0*e) + f;
        }
      }
    }
  }
//
//
// Do a check for rounding errors in the distance-squared.  It appears that
// the conventional methods for calculating sqrDist breaks down when very near
// to or at the surface (as required by transport).  We'll therefore also use
// the magnitude-squared of the vector displacement.  (Note that I've also
// tried to get around this problem by using the existing equations for
//
//    sqrDist = function(a,b,c,d,s,t)
//
// and use a more accurate addition process which minimises errors and
// breakdown of cummutitivity [where (A+B)+C != A+(B+C)] but this still
// doesn't work.
// Calculation from u = D + s*E[0] + t*E[1] is less efficient, but appears
// more robust.
//
  if (sqrDist < 0.0) { sqrDist = 0.0; }
  G4ThreeVector u = D + s*E[0] + t*E[1];
  G4double u2     = u.mag2();
//
//
// The following (part of the roundoff error check) is from Oliver Merle's
// updates.
//
  if ( sqrDist > u2 ) sqrDist = u2;

  return u;
}

///////////////////////////////////////////////////////////////////////////////
//
// Distance (G4ThreeVector, G4double)
//
// Determines the closest distance between point p and the facet.  This makes
// use of G4ThreeVector G4TriangularFacet::Distance, which stores the
// square of the distance in variable sqrDist.  If approximate methods show 
// the distance is to be greater than minDist, then forget about further
// computation and return a very large number.
//
G4double G4TriangularFacet::Distance (const G4ThreeVector &p,
  const G4double minDist)
{
//
//
// Start with quicky test to determine if the surface of the sphere enclosing
// the triangle is any closer to p than minDist.  If not, then don't bother
// about more accurate test.
//
  G4double dist = kInfinity;
  if ((p-circumcentre).mag()-radius < minDist)
  {
//
//
// It's possible that the triangle is closer than minDist, so do more accurate
// assessment.
//
    dist = Distance(p).mag();
//    dist = std::sqrt(sqrDist);
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
// forget about further computation and return a very large number.
//
// This method has been heavily modified thanks to the valuable comments and 
// corrections of Rickard Holmberg.
//
G4double G4TriangularFacet::Distance (const G4ThreeVector &p,
  const G4double minDist, const G4bool outgoing)
{
//
//
// Start with quicky test to determine if the surface of the sphere enclosing
// the triangle is any closer to p than minDist.  If not, then don't bother
// about more accurate test.
//
  G4double dist = kInfinity;
  if ((p-circumcentre).mag()-radius < minDist)
  {
//
//
// It's possible that the triangle is closer than minDist, so do more accurate
// assessment.
//
    G4ThreeVector v  = Distance(p);
    G4double dist1   = std::sqrt(sqrDist);
    G4double dir     = v.dot(surfaceNormal);
    G4bool wrongSide = (dir > 0.0 && !outgoing) || (dir < 0.0 && outgoing);
    if (dist1 <= kCarTolerance)
    {
//
//
// Point p is very close to triangle.  Check if it's on the wrong side, in
// which case return distance of 0.0 otherwise .
//
      if (wrongSide) dist = 0.0;
      else           dist = dist1;
    }
    else if (!wrongSide) dist = dist1;
  }

  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
// Extent
//
// Calculates the furthest the triangle extends in a particular direction
// defined by the vector axis.
//
G4double G4TriangularFacet::Extent (const G4ThreeVector axis)
{
  G4double s  = P0.dot(axis);
  G4double sp = P[0].dot(axis);
  if (sp > s) s = sp;
  sp = P[1].dot(axis);
  if (sp > s) s = sp;

  return s;
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
// This is to detect kSurface without doing a full Inside(p) in
// G4TessellatedSolid::Distance(p,v) calculation.
//
// This member function is thanks the valuable work of Rickard Holmberg.  PT.
// However, "gotos" are the Work of the Devil have been exorcised with
// extreme prejudice!!
//
// IMPORTANT NOTE:  These calculations are predicated on v being a unit
// vector.  If G4TessellatedSolid or other classes call this member function
// with |v| != 1 then there will be errors.
//
G4bool G4TriangularFacet::Intersect (const G4ThreeVector &p,
                   const G4ThreeVector &v, G4bool outgoing, G4double &distance,
                         G4double &distFromSurface, G4ThreeVector &normal)
{
//
//
// Check whether the direction of the facet is consistent with the vector v
// and the need to be outgoing or ingoing.  If inconsistent, disregard and
// return false.
//
  G4double w = v.dot(surfaceNormal);
  if ((outgoing && (w <-dirTolerance)) || (!outgoing && (w > dirTolerance)))
  {
    distance        = kInfinity;
    distFromSurface = kInfinity;
    normal          = G4ThreeVector(0.0,0.0,0.0);
    return false;
  }
//
//
// Calculate the orthogonal distance from p to the surface containing the
// triangle.  Then determine if we're on the right or wrong side of the
// surface (at a distance greater than kCarTolerance) to be consistent with
// "outgoing".
//
  G4ThreeVector D  = P0 - p;
  distFromSurface  = D.dot(surfaceNormal);
  G4bool wrongSide = (outgoing && (distFromSurface < -0.5*kCarTolerance)) ||
                    (!outgoing && (distFromSurface >  0.5*kCarTolerance));
  if (wrongSide)
  {
    distance        = kInfinity;
    distFromSurface = kInfinity;
    normal          = G4ThreeVector(0.0,0.0,0.0);
    return false;
  }

  wrongSide = (outgoing && (distFromSurface < 0.0)) ||
             (!outgoing && (distFromSurface > 0.0));
  if (wrongSide)
  {
//
//
// We're slightly on the wrong side of the surface.  Check if we're close
// enough using a precise distance calculation.
//
    G4ThreeVector u = Distance(p);
    if ( sqrDist <= kCarTolerance*kCarTolerance )
    {
//
//
// We're very close.  Therefore return a small negative number to pretend
// we intersect.
//
//      distance = -0.5*kCarTolerance;
      distance = 0.0;
      normal   = surfaceNormal;
      return true;
    }
    else
    {
//
//
// We're close to the surface containing the triangle, but sufficiently
// far from the triangle, and on the wrong side compared to the directions
// of the surface normal and v.  There is no intersection.
//
      distance        = kInfinity;
      distFromSurface = kInfinity;
      normal          = G4ThreeVector(0.0,0.0,0.0);
      return false;
    }
  }
  if (w < dirTolerance && w > -dirTolerance)
  {
//
//
// The ray is within the plane of the triangle.  Project the problem into 2D
// in the plane of the triangle.  First try to create orthogonal unit vectors
// mu and nu, where mu is E[0]/|E[0]|.  This is kinda like
// the original algorithm due to Rickard Holmberg, but with better mathematical
// justification than the original method ... however, beware Rickard's was less
// time-consuming.
//
// Note that vprime is not a unit vector.  We need to keep it unnormalised
// since the values of distance along vprime (s0 and s1) for intersection with
// the triangle will be used to determine if we cut the plane at the same
// time.
//
    G4ThreeVector mu = E[0].unit();
    G4ThreeVector nu = surfaceNormal.cross(mu);
    G4TwoVector pprime(p.dot(mu),p.dot(nu));
    G4TwoVector vprime(v.dot(mu),v.dot(nu));
    G4TwoVector P0prime(P0.dot(mu),P0.dot(nu));
    G4TwoVector E0prime(E[0].mag(),0.0);
    G4TwoVector E1prime(E[1].dot(mu),E[1].dot(nu));

    G4TwoVector loc[2];
    if ( tGeomAlg->IntersectLineAndTriangle2D(pprime,vprime,P0prime,
                                              E0prime,E1prime,loc) )
    {
//
//
// There is an intersection between the line and triangle in 2D.  Now check
// which part of the line intersects with the plane containing the triangle
// in 3D.
//
      G4double vprimemag = vprime.mag();
      G4double s0        = (loc[0] - pprime).mag()/vprimemag;
      G4double s1        = (loc[1] - pprime).mag()/vprimemag;
      G4double normDist0 = surfaceNormal.dot(s0*v) - distFromSurface;
      G4double normDist1 = surfaceNormal.dot(s1*v) - distFromSurface;

      if ((normDist0 < 0.0 && normDist1 < 0.0) ||
          (normDist0 > 0.0 && normDist1 > 0.0) ||
          (normDist0 == 0.0 && normDist1 == 0.0) ) 
      {
        distance        = kInfinity;
        distFromSurface = kInfinity;
        normal          = G4ThreeVector(0.0,0.0,0.0);
        return false;
      }
      else
      {
        G4double dnormDist = normDist1-normDist0;
        if (std::fabs(dnormDist) < DBL_EPSILON)
        {
          distance = s0;
          normal   = surfaceNormal;
          if (!outgoing) distFromSurface = -distFromSurface;
          return true;
        }
        else
        {
          distance = s0 - normDist0*(s1-s0)/dnormDist;
          normal   = surfaceNormal;
          if (!outgoing) distFromSurface = -distFromSurface;
          return true;
        }
      }

//      G4ThreeVector dloc   = loc1 - loc0;
//      G4ThreeVector dlocXv = dloc.cross(v);
//      G4double dlocXvmag   = dlocXv.mag();
//      if (dloc.mag() <= 0.5*kCarTolerance || dlocXvmag <= DBL_EPSILON)
//      {
//        distance = loc0.mag();
//        normal = surfaceNormal;
//        if (!outgoing) distFromSurface = -distFromSurface;
//        return true;
//      }

//      G4ThreeVector loc0Xv   = loc0.cross(v);
//      G4ThreeVector loc1Xv   = loc1.cross(v);
//      G4double sameDir       = -loc0Xv.dot(loc1Xv);
//      if (sameDir < 0.0)
//      {
//        distance        = kInfinity;
//        distFromSurface = kInfinity;
//        normal          = G4ThreeVector(0.0,0.0,0.0);
//        return false;
//      }
//      else
//      {
//        distance = loc0.mag() + loc0Xv.mag() * dloc.mag()/dlocXvmag;
//        normal   = surfaceNormal;
//        if (!outgoing) distFromSurface = -distFromSurface;
//        return true;
//      }
    }
    else
    {
      distance        = kInfinity;
      distFromSurface = kInfinity;
      normal          = G4ThreeVector(0.0,0.0,0.0);
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
  distance         = distFromSurface / w;
  G4ThreeVector pp = p + v*distance;
  G4ThreeVector DD = P0 - pp;
  G4double d       = E[0].dot(DD);
  G4double e       = E[1].dot(DD);
  G4double s       = b*e - c*d;
  G4double t       = b*d - a*e;

  G4double sTolerance = (std::fabs(b)+ std::fabs(c) + std::fabs(d)
                       + std::fabs(e)) *kCarTolerance;
  G4double tTolerance = (std::fabs(a)+ std::fabs(b) + std::fabs(d)
                       + std::fabs(e)) *kCarTolerance;
  G4double detTolerance = (std::fabs(a)+ std::fabs(c)
                       + 2*std::fabs(b) ) *kCarTolerance;

  //if (s < 0.0 || t < 0.0 || s+t > det)
  if (s < -sTolerance || t < -tTolerance || ( s+t - det ) > detTolerance)
  {
//
//
// The intersection is outside of the triangle.
//
    distance        = kInfinity;
    distFromSurface = kInfinity;
    normal          = G4ThreeVector(0.0,0.0,0.0);
    return false;
  }
  else
  {
//
//
// There is an intersection.  Now we only need to set the surface normal.
//
     normal = surfaceNormal;
     if (!outgoing) distFromSurface = -distFromSurface;
     return true;
  }
}

////////////////////////////////////////////////////////////////////////
//
// GetPointOnFace
//
// Auxiliary method for get a random point on surface

G4ThreeVector G4TriangularFacet::GetPointOnFace() const
{
  G4double alpha = G4RandFlat::shoot(0.,1.);
  G4double beta = G4RandFlat::shoot(0.,1);
  G4double lambda1=alpha*beta;
  G4double lambda0=alpha-lambda1;
  
  return (P0 + lambda0*E[0] + lambda1*E[1]);
}

////////////////////////////////////////////////////////////////////////
//
// GetArea
//
// Auxiliary method for returning the surface area

G4double G4TriangularFacet::GetArea()
{
  return area;
}
////////////////////////////////////////////////////////////////////////
//
