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
//
// $Id: G4QuadrangularFacet.cc,v 1.9 2010-09-23 10:27:25 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4QuadrangularFacet.cc
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

#include "G4QuadrangularFacet.hh"
#include "globals.hh"
#include "Randomize.hh"

///////////////////////////////////////////////////////////////////////////////
//
// !!!THIS IS A FUDGE!!!  IT'S TWO ADJACENT G4TRIANGULARFACETS
// --- NOT EFFICIENT BUT PRACTICAL.
//
G4QuadrangularFacet::G4QuadrangularFacet (const G4ThreeVector Pt0,
                 const G4ThreeVector vt1, const G4ThreeVector vt2,
                 const G4ThreeVector vt3, G4FacetVertexType vertexType)
  : G4VFacet(), facet1(0), facet2(0)
{
  P0        = Pt0;
  nVertices = 4;
  if (vertexType == ABSOLUTE)
  {
    P.push_back(vt1);
    P.push_back(vt2);
    P.push_back(vt3);
  
    E.push_back(vt1 - P0);
    E.push_back(vt2 - P0);
    E.push_back(vt3 - P0);
  }
  else
  {
    P.push_back(P0 + vt1);
    P.push_back(P0 + vt2);
    P.push_back(P0 + vt3);
  
    E.push_back(vt1);
    E.push_back(vt2);
    E.push_back(vt3);
  }

  G4double length1 = E[0].mag();
  G4double length2 = (P[1]-P[0]).mag();
  G4double length3 = (P[2]-P[1]).mag();
  G4double length4 = E[2].mag();
  
  G4ThreeVector normal1 = E[0].cross(E[1]).unit();
  G4ThreeVector normal2 = E[1].cross(E[2]).unit(); 
  
  if (length1 <= kCarTolerance || length2 <= kCarTolerance ||
      length3 <= kCarTolerance || length4 <= kCarTolerance ||
      normal1.dot(normal2) < 0.9999999999)
  {
    G4Exception("G4QuadrangularFacet::G4QuadrangularFacet()",
                "InvalidSetup", JustWarning,
                "Length of sides of facet are too small or sides not planar.");
    G4cerr << G4endl;
    G4cerr << "P0 = " << P0   << G4endl;
    G4cerr << "P1 = " << P[0] << G4endl;
    G4cerr << "P2 = " << P[1] << G4endl;
    G4cerr << "P3 = " << P[2] << G4endl;
    G4cerr << "Side lengths = P0->P1" << length1 << G4endl;    
    G4cerr << "Side lengths = P1->P2" << length2 << G4endl;    
    G4cerr << "Side lengths = P2->P3" << length3 << G4endl;    
    G4cerr << "Side lengths = P3->P0" << length4 << G4endl;    
    G4cerr << G4endl;
    
    isDefined     = false;
    geometryType  = "G4QuadragularFacet";
    surfaceNormal = G4ThreeVector(0.0,0.0,0.0);
  }
  else
  {
    isDefined     = true;
    geometryType  = "G4QuadrangularFacet";
    
    facet1 = new G4TriangularFacet(P0,P[0],P[1],ABSOLUTE);
    facet2 = new G4TriangularFacet(P0,P[1],P[2],ABSOLUTE);
    surfaceNormal = normal1;
    
    G4ThreeVector vtmp = 0.5 * (E[0] + E[1]);
    circumcentre       = P0 + vtmp;
    radiusSqr          = vtmp.mag2();
    radius             = std::sqrt(radiusSqr);
  
    for (size_t i=0; i<4; i++) I.push_back(0);
  }
}

///////////////////////////////////////////////////////////////////////////////
//
G4QuadrangularFacet::~G4QuadrangularFacet ()
{
  delete facet1;
  delete facet2;
  
  P.clear();
  E.clear();
  I.clear();
}

///////////////////////////////////////////////////////////////////////////////
//
G4QuadrangularFacet::G4QuadrangularFacet (const G4QuadrangularFacet &rhs)
  : G4VFacet(rhs)
{
  facet1 = new G4TriangularFacet(*(rhs.facet1));
  facet2 = new G4TriangularFacet(*(rhs.facet2));
}

///////////////////////////////////////////////////////////////////////////////
//
const G4QuadrangularFacet &
G4QuadrangularFacet::operator=(G4QuadrangularFacet &rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VFacet::operator=(rhs);

   // Copy data
   //
   delete facet1; facet1 = new G4TriangularFacet(*(rhs.facet1));
   delete facet2; facet2 = new G4TriangularFacet(*(rhs.facet2));

   return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
G4VFacet *G4QuadrangularFacet::GetClone ()
{
  G4QuadrangularFacet *c =
    new G4QuadrangularFacet (P0, P[0], P[1], P[2], ABSOLUTE);
  G4VFacet *cc         = 0;
  cc                   = c;
  return cc;
}

///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector G4QuadrangularFacet::Distance (const G4ThreeVector &p)
{
  G4ThreeVector v1 = facet1->Distance(p);
  G4ThreeVector v2 = facet2->Distance(p);
  
  if (v1.mag2() < v2.mag2()) return v1;
  else return v2;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4QuadrangularFacet::Distance (const G4ThreeVector &p,
  const G4double)
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
  
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4QuadrangularFacet::Distance (const G4ThreeVector &p,
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
G4double G4QuadrangularFacet::Extent (const G4ThreeVector axis)
{
  G4double s  = P0.dot(axis);
  for (G4ThreeVectorList::iterator it=P.begin(); it!=P.end(); it++)
  {
    G4double sp = it->dot(axis);
    if (sp > s) s = sp;
  }

  return s;
}

///////////////////////////////////////////////////////////////////////////////
//
G4bool G4QuadrangularFacet::Intersect (const G4ThreeVector &p,
  const G4ThreeVector &v, G4bool outgoing, G4double &distance,
  G4double &distFromSurface, G4ThreeVector &normal)
{
  G4bool intersect =
    facet1->Intersect(p,v,outgoing,distance,distFromSurface,normal);
  if (!intersect)
  {
    intersect = facet2->Intersect(p,v,outgoing,distance,distFromSurface,normal);
  }
  
  if (!intersect)
  {
    distance        = kInfinity;
    distFromSurface = kInfinity;
    normal          = G4ThreeVector(0.0,0.0,0.0);
  }
  
  return intersect;
}

////////////////////////////////////////////////////////////////////////
//
// GetPointOnFace
//
// Auxiliary method for get a random point on surface

G4ThreeVector G4QuadrangularFacet::GetPointOnFace() const
{
  G4ThreeVector pr;

  if ( G4UniformRand() < 0.5 )
  {
    pr = facet1->GetPointOnFace();
  }
  else
  {
    pr = facet2->GetPointOnFace();
  }

  return pr;
}

////////////////////////////////////////////////////////////////////////
//
// GetArea
//
// Auxiliary method for returning the surface area

G4double G4QuadrangularFacet::GetArea()
{
  if (!area)  { area = facet1->GetArea() + facet2->GetArea(); }

  return area;
}
