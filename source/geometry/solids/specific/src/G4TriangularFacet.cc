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
// * subject DEFCON 705 IPR conditions.                               *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4TriangularFacet.cc,v 1.5 2006/06/29 18:49:02 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "G4TriangularFacet.hh"
#include "globals.hh"

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
  : G4VFacet()
{
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
    G4Exception("G4TriangularFacet::G4TriangularFacet()", "InvalidSetup",
                JustWarning, "Length of sides of facet are too small.");
    G4cerr << G4endl;
    G4cerr << "P0 = " << P0   << G4endl;
    G4cerr << "P1 = " << P[0] << G4endl;
    G4cerr << "P2 = " << P[1] << G4endl;
    G4cerr << "Side lengths = P0->P1" << Emag1 << G4endl;    
    G4cerr << "Side lengths = P0->P2" << Emag2 << G4endl;    
    G4cerr << "Side lengths = P1->P2" << Emag3 << G4endl;    
    G4cerr << G4endl;
    
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
    det = std::abs(a*c - b*b);
    
    sMin = -0.5*kCarTolerance/std::sqrt(a);
    sMax = 1.0 - sMin;
    tMin = -0.5*kCarTolerance/std::sqrt(c);
    
    G4ThreeVector vtmp = 0.25 * (E[0] + E[1]);
    centroid           = P0 + vtmp;
    radiusSqr          = vtmp.mag2();
    radius             = std::sqrt(radiusSqr);
  
    for (size_t i=0; i<3; i++) I.push_back(0);
  }
}

///////////////////////////////////////////////////////////////////////////////
//
G4TriangularFacet::~G4TriangularFacet ()
{
  P.clear();
  E.clear();
  I.clear();
}

///////////////////////////////////////////////////////////////////////////////
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
G4TriangularFacet *G4TriangularFacet::GetFlippedFacet ()
{
  G4TriangularFacet *flipped = new G4TriangularFacet (P0, P[1], P[0], ABSOLUTE);
  return flipped;
}

///////////////////////////////////////////////////////////////////////////////
//
// Determine the closest distance from the facet to the point p.  If the
// direction of the vector to the closest point is outward-going and outgoing
// is true or the vector is in-going and outgoing is false then the distance
// is returned.  Otherwise kInfinity is returned.
//
G4ThreeVector G4TriangularFacet::Distance (const G4ThreeVector &p)
{
  G4ThreeVector D  = P0 - p;
  G4double d       = E[0].dot(D);
  G4double e       = E[1].dot(D);
  G4double f       = D.mag2();
  G4double s       = b*e - c*d;
  G4double t       = b*d - a*e;
  G4double sqrDist = 0.0;
  
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
        else if (-e >= c)  {t = 0.0; sqrDist = c + 2.0*e + f;}
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
      G4double invDet = 1.0 / det;
      s              *= invDet;
      t              *= invDet;
      sqrDist         = s*(a*s + b*t + 2.0*d) + t*(b*s + c*t + 2.0*e) + f;
    }
  }
  else
  {
    G4double tmp0  = 0.0;
    G4double tmp1  = 0.0;
    G4double numer = 0.0;
    G4double denom = 0.0;
    if (s < 0.0)
    {
   //
   // We are in region 2.
   //
      tmp0 = b + d;
      tmp1 = c + e;
      if (tmp1 > tmp0)
      {
        numer = tmp1 - tmp0;
        denom = a - 2.0*b*c;
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
      tmp0 = b + e;
      tmp1 = a + d;
      if (tmp1 > tmp0)
      {
        numer = tmp1 - tmp0;
        denom = a - 2.0*b*c;
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
      numer = c + f - b - d;
      if (numer <= 0.0)
      {
        s       = 0.0;
        t       = 1.0;
        sqrDist = c + 2.0*e*f;
      }
      else
      {
        denom = a - 2.0*b*c;
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
  return D + s*E[0] + t*E[1];
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TriangularFacet::Distance (const G4ThreeVector &p,
  const G4double minDist)
{
  /*G4ThreeVector D  = P0 - p;
  G4double d       = E[0].dot(D);
  G4double e       = E[1].dot(D);
  G4double s       = b*e - c*d;
  G4double t       = b*d - a*e;*/
  G4double dist = kInfinity;
  
  /*if (s+t > 1.0 || s < 0.0 || t < 0.0)
  {
    G4ThreeVector D0 = P0 - p;
    G4ThreeVector D1 = P[0] - p;
    G4ThreeVector D2 = P[1] - p;
    
    G4double d0 = D0.mag();
    G4double d1 = D1.mag();
    G4double d2 = D2.mag();
    
    dist = min(d0, min(d1, d2));
    if (dist > minDist) return kInfinity;
  }*/
  
  dist = Distance(p).mag();
  if (dist > minDist) return kInfinity;
  
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
// Determine the distance to point p bearing in mind that if the distance is
// likely to be longer than minDist, forget doing further calculation and
// return kInfinity.
//
G4double G4TriangularFacet::Distance (const G4ThreeVector &p,
                      const G4double, const G4bool outgoing)
{
  /*G4ThreeVector D  = P0 - p;
  G4double d       = E[0].dot(D);
  G4double e       = E[1].dot(D);
  G4double s       = b*e - c*d;
  G4double t       = b*d - a*e;*/
  G4double dist = kInfinity;
  
  /*if (s+t > 1.0 || s < 0.0 || t < 0.0)
  {
    G4ThreeVector D0 = P0 - p;
    G4ThreeVector D1 = P[0] - p;
    G4ThreeVector D2 = P[1] - p;
    
    G4double d0 = D0.mag();
    G4double d1 = D1.mag();
    G4double d2 = D2.mag();
    
    dist = min(d0, min(d1, d2));
    if (dist > minDist ||
      (D0.dot(surfaceNormal) > 0.0 && !outgoing) ||
      (D0.dot(surfaceNormal) < 0.0 && outgoing)) return kInfinity;
  }*/
  
  G4ThreeVector v = Distance(p);
  G4double dir    = v.dot(surfaceNormal);
  if ((dir > dirTolerance && !outgoing) ||
      (dir <-dirTolerance && outgoing)) dist = kInfinity;
  else dist = v.mag();
  
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
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
G4bool G4TriangularFacet::Intersect (const G4ThreeVector &p,
                   const G4ThreeVector &v, G4bool outgoing, G4double &distance,
                         G4double &distFromSurface, G4ThreeVector &normal)
{
  G4ThreeVector D = P0 - p;
  G4double d      = E[0].dot(D);
  G4double e      = E[1].dot(D);
  G4double g      = E[0].dot(v);
  G4double h      = E[1].dot(v);
  G4double q      = D.dot(v);
  
  G4double A00  = a - g*g;
  G4double A11  = c - h*h;
  G4double A01  = b - g*h;
  G4double det2 = A00*A11 - A01*A01;
  
  G4double s = kInfinity;
  G4double t = kInfinity;
  
  G4double dist       = kInfinity;
  G4bool intersect    = false;
  G4double normalComp = 0.0;
  
  
  if (det2 != 0.0)
  {
    G4double B0 = q*g - d;
    G4double B1 = q*h - e;
    s           = (A11*B0 - A01*B1)/det2;
    if ((s >= sMin) && (s <= sMax))
    {
      t = (A00*B1 - A01*B0)/det2;
      if ((t >= tMin) && (t < 1.0 - s + std::fabs(sMin)))
      {                                //THIS IS A FUDGE FOR THE MOMENT
        dist       = q + g*s + h*t;
        normalComp = v.dot(surfaceNormal);
//      intersect  = (dist >= 0.0 &&
//        ((outgoing && normalComp > 0.0) || (!outgoing && normalComp < 0.0)));
        intersect  = (dist >= -kCarTolerance*0.5 &&
          ((outgoing && normalComp > dirTolerance) ||
          (!outgoing && normalComp <-dirTolerance))); //FUDGE FOR THE MOMENT
      }
    }
  }
  
  if (intersect)
  {
    if (dist > kCarTolerance * 0.5) distance = dist;
    else dist = 0.0;
    distFromSurface = dist * normalComp;
    normal          = surfaceNormal;
  }
  else
  {
    distance        = kInfinity;
    distFromSurface = kInfinity;
    normal          = G4ThreeVector(0.0,0.0,0.0);
  }
  
  return intersect;
}
