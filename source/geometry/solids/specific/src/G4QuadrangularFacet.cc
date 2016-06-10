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
// $Id: G4QuadrangularFacet.cc 67011 2013-01-29 16:17:41Z gcosmo $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK - Created.
// 12 October 2012, M Gayer, CERN
//                  New implementation reducing memory requirements by 50%,
//                  and considerable CPU speedup together with the new
//                  implementation of G4TessellatedSolid.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "G4QuadrangularFacet.hh"
#include "geomdefs.hh"
#include "Randomize.hh"
 
using namespace std;

///////////////////////////////////////////////////////////////////////////////
//
// !!!THIS IS A FUDGE!!!  IT'S TWO ADJACENT G4TRIANGULARFACETS
// --- NOT EFFICIENT BUT PRACTICAL.
//
G4QuadrangularFacet::G4QuadrangularFacet (const G4ThreeVector &vt0,
                                          const G4ThreeVector &vt1,
                                          const G4ThreeVector &vt2,
                                          const G4ThreeVector &vt3,
                                                G4FacetVertexType vertexType)
{
  G4ThreeVector e1, e2, e3;

  SetVertex(0, vt0);
  if (vertexType == ABSOLUTE)
  {
    SetVertex(1, vt1);
    SetVertex(2, vt2);
    SetVertex(3, vt3);

    e1 = vt1 - vt0;
    e2 = vt2 - vt0;
    e3 = vt3 - vt0;
  }
  else
  {
    SetVertex(1, vt0 + vt1);
    SetVertex(2, vt0 + vt2);
    SetVertex(3, vt0 + vt3);

    e1 = vt1;
    e2 = vt2;
    e3 = vt3;
  }
  G4double length1 = e1.mag();
  G4double length2 = (GetVertex(2)-GetVertex(1)).mag();
  G4double length3 = (GetVertex(3)-GetVertex(2)).mag();
  G4double length4 = e3.mag();

  G4ThreeVector normal1 = e1.cross(e2).unit();
  G4ThreeVector normal2 = e2.cross(e3).unit(); 

  bool isDefined = (length1 > kCarTolerance && length2 > kCarTolerance &&
    length3 > kCarTolerance && length4 > kCarTolerance &&
    normal1.dot(normal2) >= 0.9999999999);

  if (isDefined)
  {
    fFacet1 = G4TriangularFacet (GetVertex(0),GetVertex(1),
                                 GetVertex(2),ABSOLUTE);
    fFacet2 = G4TriangularFacet (GetVertex(0),GetVertex(2),
                                 GetVertex(3),ABSOLUTE);

    G4TriangularFacet facet3 (GetVertex(0),GetVertex(1),GetVertex(3),ABSOLUTE);
    G4TriangularFacet facet4 (GetVertex(1),GetVertex(2),GetVertex(3),ABSOLUTE);

    G4ThreeVector normal12 = fFacet1.GetSurfaceNormal()
                           + fFacet2.GetSurfaceNormal();
    G4ThreeVector normal34 = facet3.GetSurfaceNormal()
                           + facet4.GetSurfaceNormal();
    G4ThreeVector normal = 0.25 * (normal12 + normal34);

    fFacet1.SetSurfaceNormal (normal);
    fFacet2.SetSurfaceNormal (normal);

    G4ThreeVector vtmp = 0.5 * (e1 + e2);
    fCircumcentre = GetVertex(0) + vtmp;
    G4double radiusSqr = vtmp.mag2();
    fRadius = std::sqrt(radiusSqr);
  }
  else
  {
    G4Exception("G4QuadrangularFacet::G4QuadrangularFacet()",
                "GeomSolids0002", JustWarning,
                "Length of sides of facet are too small or sides not planar.");
    G4cout << G4endl;
    G4cout << "P0 = " << GetVertex(0) << G4endl;
    G4cout << "P1 = " << GetVertex(1) << G4endl;
    G4cout << "P2 = " << GetVertex(2) << G4endl;
    G4cout << "P3 = " << GetVertex(3) << G4endl;
    G4cout << "Side lengths = P0->P1" << length1 << G4endl;    
    G4cout << "Side lengths = P1->P2" << length2 << G4endl;    
    G4cout << "Side lengths = P2->P3" << length3 << G4endl;    
    G4cout << "Side lengths = P3->P0" << length4 << G4endl;    
    G4cout << G4endl;
    fRadius = 0.0;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
G4QuadrangularFacet::~G4QuadrangularFacet ()
{
}

///////////////////////////////////////////////////////////////////////////////
//
G4QuadrangularFacet::G4QuadrangularFacet (const G4QuadrangularFacet &rhs)
  : G4VFacet(rhs)
{
  fFacet1 = rhs.fFacet1;
  fFacet2 = rhs.fFacet2;
  fRadius = 0.0;
}

///////////////////////////////////////////////////////////////////////////////
//
G4QuadrangularFacet &
G4QuadrangularFacet::operator=(const G4QuadrangularFacet &rhs)
{
  if (this == &rhs)
    return *this;

  fFacet1 = rhs.fFacet1;
  fFacet2 = rhs.fFacet2;
  fRadius = 0.0;

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
G4VFacet *G4QuadrangularFacet::GetClone ()
{
  G4QuadrangularFacet *c = new G4QuadrangularFacet (GetVertex(0), GetVertex(1),
                                                    GetVertex(2), GetVertex(3),
                                                    ABSOLUTE);
  return c;
}

///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector G4QuadrangularFacet::Distance (const G4ThreeVector &p)
{
  G4ThreeVector v1 = fFacet1.Distance(p);
  G4ThreeVector v2 = fFacet2.Distance(p);

  if (v1.mag2() < v2.mag2()) return v1;
  else return v2;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4QuadrangularFacet::Distance (const G4ThreeVector &p,
                                              G4double)
{  
  G4double dist = Distance(p).mag();
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4QuadrangularFacet::Distance (const G4ThreeVector &p, G4double,
                                        const G4bool outgoing)
{
  G4double dist;

  G4ThreeVector v = Distance(p);
  G4double dir = v.dot(GetSurfaceNormal());
  if ( ((dir > dirTolerance) && (!outgoing))
    || ((dir < -dirTolerance) && outgoing))
    dist = kInfinity;
  else 
    dist = v.mag();
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4QuadrangularFacet::Extent (const G4ThreeVector axis)
{
  G4double ss  = 0;

  for (G4int i = 0; i <= 3; ++i)
  {
    G4double sp = GetVertex(i).dot(axis);
    if (sp > ss) ss = sp;
  }
  return ss;
}

///////////////////////////////////////////////////////////////////////////////
//
G4bool G4QuadrangularFacet::Intersect (const G4ThreeVector &p,
                                       const G4ThreeVector &v,
                                             G4bool outgoing,
                                             G4double &distance,
                                             G4double &distFromSurface,
                                             G4ThreeVector &normal)
{
  G4bool intersect =
    fFacet1.Intersect(p,v,outgoing,distance,distFromSurface,normal);
  if (!intersect) intersect =
    fFacet2.Intersect(p,v,outgoing,distance,distFromSurface,normal);
  if (!intersect)
  {
    distance = distFromSurface = kInfinity;
    normal.set(0,0,0);
  }
  return intersect;
}

///////////////////////////////////////////////////////////////////////////////
//
// Auxiliary method to get a random point on surface
//
G4ThreeVector G4QuadrangularFacet::GetPointOnFace() const
{
  G4ThreeVector pr = (G4RandFlat::shoot(0.,1.) < 0.5)
                   ? fFacet1.GetPointOnFace() : fFacet2.GetPointOnFace();
  return pr;
}

///////////////////////////////////////////////////////////////////////////////
//
// Auxiliary method for returning the surface area
//
G4double G4QuadrangularFacet::GetArea()
{
  G4double area = fFacet1.GetArea() + fFacet2.GetArea();
  return area;
}

///////////////////////////////////////////////////////////////////////////////
//
G4String G4QuadrangularFacet::GetEntityType () const
{
  return "G4QuadrangularFacet";
}

///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector G4QuadrangularFacet::GetSurfaceNormal () const
{
  return fFacet1.GetSurfaceNormal();
}
