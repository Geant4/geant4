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
// G4QuadrangularFacet class implementation.
//
// 31 October 2004 P R Truscott, QinetiQ Ltd, UK - Created.
// 12 October 2012 M Gayer, CERN
//                 New implementation reducing memory requirements by 50%,
//                 and considerable CPU speedup together with the new
//                 implementation of G4TessellatedSolid.
// 29 February 2016 E Tcherniaev, CERN
//                 Added exhaustive tests to catch various problems with a
//                 quadrangular facet: collinear vertices, non planar surface,
//                 degenerate, concave or self intersecting quadrilateral.
// --------------------------------------------------------------------

#include "G4QuadrangularFacet.hh"
#include "geomdefs.hh"
#include "Randomize.hh"
 
using namespace std;

///////////////////////////////////////////////////////////////////////////////
//
// Constructing two adjacent G4TriangularFacet
// Not efficient, but practical...
//
G4QuadrangularFacet::G4QuadrangularFacet (const G4ThreeVector& vt0,
                                          const G4ThreeVector& vt1,
                                          const G4ThreeVector& vt2,
                                          const G4ThreeVector& vt3,
                                                G4FacetVertexType vertexType)
  : G4VFacet()
{
  G4double delta   =  1.0 * kCarTolerance; // dimension tolerance
  G4double epsilon = 0.01 * kCarTolerance; // planarity tolerance

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

  // Check length of sides and diagonals
  //
  G4double leng1 = e1.mag();
  G4double leng2 = (e2-e1).mag();
  G4double leng3 = (e3-e2).mag();
  G4double leng4 = e3.mag();

  G4double diag1 = e2.mag();
  G4double diag2 = (e3-e1).mag();

  if (leng1 <= delta || leng2 <= delta || leng3 <= delta || leng4 <= delta ||
      diag1 <= delta || diag2 <= delta)
  {
    ostringstream message;
    message << "Sides/diagonals of facet are too small." << G4endl
            << "P0 = " << GetVertex(0) << G4endl
            << "P1 = " << GetVertex(1) << G4endl
            << "P2 = " << GetVertex(2) << G4endl
            << "P3 = " << GetVertex(3) << G4endl
            << "Side1 length (P0->P1) = " << leng1 << G4endl
            << "Side2 length (P1->P2) = " << leng2 << G4endl
            << "Side3 length (P2->P3) = " << leng3 << G4endl
            << "Side4 length (P3->P0) = " << leng4 << G4endl
            << "Diagonal1 length (P0->P2) = " << diag1 << G4endl
            << "Diagonal2 length (P1->P3) = " << diag2;
    G4Exception("G4QuadrangularFacet::G4QuadrangularFacet()",
                "GeomSolids1001", JustWarning, message);
    return;
  }

  // Check that vertices are not collinear
  //
  G4double s1 = (e1.cross(e2)).mag()*0.5;
  G4double s2 = ((e2-e1).cross(e3-e2)).mag()*0.5;
  G4double s3 = (e2.cross(e3)).mag()*0.5;
  G4double s4 = (e1.cross(e3)).mag()*0.5;

  G4double h1 = 2.*s1 / std::max(std::max(leng1,leng2),diag1);
  G4double h2 = 2.*s2 / std::max(std::max(leng2,leng3),diag2);
  G4double h3 = 2.*s3 / std::max(std::max(leng3,leng4),diag1);
  G4double h4 = 2.*s4 / std::max(std::max(leng4,leng1),diag2);
   
  if (h1 <= delta || h2 <= delta || h3 <= delta || h4 <= delta )
  {
    ostringstream message;
    message << "Facet has three or more collinear vertices." << G4endl
            << "P0 = " << GetVertex(0) << G4endl
            << "P1 = " << GetVertex(1) << G4endl
            << "P2 = " << GetVertex(2) << G4endl
            << "P3 = " << GetVertex(3) << G4endl
            << "Smallest heights:" << G4endl
            << "  in triangle P0-P1-P2 = " << h1 << G4endl
            << "  in triangle P1-P2-P3 = " << h2 << G4endl
            << "  in triangle P2-P3-P0 = " << h3 << G4endl
            << "  in triangle P3-P0-P1 = " << h4;
    G4Exception("G4QuadrangularFacet::G4QuadrangularFacet()",
	        "GeomSolids1001", JustWarning, message);
    return;
  }

  // Check that vertices are coplanar by computing minimal
  // height of tetrahedron comprising of vertices
  //
  G4double smax = std::max( std::max(s1,s2), std::max(s3,s4) ); 
  G4double hmin = 0.5 * std::fabs( e1.dot(e2.cross(e3)) ) / smax;
  if (hmin >= epsilon)
  {
    ostringstream message;
    message << "Facet is not planar." << G4endl
            << "Disrepancy = " << hmin << G4endl
            << "P0 = " << GetVertex(0) << G4endl
            << "P1 = " << GetVertex(1) << G4endl
            << "P2 = " << GetVertex(2) << G4endl
            << "P3 = " << GetVertex(3);
    G4Exception("G4QuadrangularFacet::G4QuadrangularFacet()",
	        "GeomSolids1001", JustWarning, message);
    return;
  }
  
  // Check that facet is convex by computing crosspoint 
  // of diagonals
  //
  G4ThreeVector normal = e2.cross(e3-e1);
  G4double s = kInfinity, t = kInfinity, magnitude2 = normal.mag2();
  if (magnitude2 > delta*delta) // check: magnitude2 != 0.
  {
    s = normal.dot(e1.cross(e3-e1)) / magnitude2;
    t = normal.dot(e1.cross(e2))    / magnitude2;
  }
  if (s <= 0. || s >= 1. || t <= 0. || t >= 1.)
  {
    ostringstream message;
    message << "Facet is not convex." << G4endl
            << "Parameters of crosspoint of diagonals: "
            << s << " and " << t << G4endl
            << "should both be within (0,1) range" << G4endl
	    << "P0 = " << GetVertex(0) << G4endl
	    << "P1 = " << GetVertex(1) << G4endl
	    << "P2 = " << GetVertex(2) << G4endl
	    << "P3 = " << GetVertex(3);
    G4Exception("G4QuadrangularFacet::G4QuadrangularFacet()",
	        "GeomSolids1001", JustWarning, message);
    return;
  }

  // Define facet
  //
  fFacet1 = G4TriangularFacet(GetVertex(0),GetVertex(1),GetVertex(2),ABSOLUTE);
  fFacet2 = G4TriangularFacet(GetVertex(0),GetVertex(2),GetVertex(3),ABSOLUTE);

  normal = normal.unit();
  fFacet1.SetSurfaceNormal(normal);
  fFacet2.SetSurfaceNormal(normal);
  
  G4ThreeVector vtmp = 0.5 * (e1 + e2);
  fCircumcentre = GetVertex(0) + vtmp;
  G4double radiusSqr = vtmp.mag2();
  fRadius = std::sqrt(radiusSqr);
  // 29.02.2016 Remark by E.Tcherniaev: computation
  // of fCircumcenter and fRadius is wrong, however
  // it did not create any problem till now.
  // Bizarre! Need to investigate!
}

///////////////////////////////////////////////////////////////////////////////
//
G4QuadrangularFacet::~G4QuadrangularFacet ()
{
}

///////////////////////////////////////////////////////////////////////////////
//
G4QuadrangularFacet::G4QuadrangularFacet (const G4QuadrangularFacet& rhs)
  : G4VFacet(rhs)
{
  fFacet1 = rhs.fFacet1;
  fFacet2 = rhs.fFacet2;
  fRadius = 0.0;
}

///////////////////////////////////////////////////////////////////////////////
//
G4QuadrangularFacet &
G4QuadrangularFacet::operator=(const G4QuadrangularFacet& rhs)
{
  if (this == &rhs)  return *this;

  fFacet1 = rhs.fFacet1;
  fFacet2 = rhs.fFacet2;
  fRadius = 0.0;

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
G4VFacet* G4QuadrangularFacet::GetClone ()
{
  G4QuadrangularFacet *c = new G4QuadrangularFacet (GetVertex(0), GetVertex(1),
                                                    GetVertex(2), GetVertex(3),
                                                    ABSOLUTE);
  return c;
}

///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector G4QuadrangularFacet::Distance (const G4ThreeVector& p)
{
  G4ThreeVector v1 = fFacet1.Distance(p);
  G4ThreeVector v2 = fFacet2.Distance(p);

  if (v1.mag2() < v2.mag2()) return v1;
  else return v2;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4QuadrangularFacet::Distance (const G4ThreeVector& p,
                                              G4double)
{  
  G4double dist = Distance(p).mag();
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4QuadrangularFacet::Distance (const G4ThreeVector& p, G4double,
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
G4bool G4QuadrangularFacet::Intersect (const G4ThreeVector& p,
                                       const G4ThreeVector& v,
                                             G4bool outgoing,
                                             G4double& distance,
                                             G4double& distFromSurface,
                                             G4ThreeVector& normal)
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
// Auxiliary method to get a uniform random point on the facet
//
G4ThreeVector G4QuadrangularFacet::GetPointOnFace() const
{
  G4double s1 = fFacet1.GetArea();
  G4double s2 = fFacet2.GetArea();
  return ((s1+s2)*G4UniformRand() < s1) ?
    fFacet1.GetPointOnFace() : fFacet2.GetPointOnFace();
}

///////////////////////////////////////////////////////////////////////////////
//
// Auxiliary method for returning the surface area
//
G4double G4QuadrangularFacet::GetArea() const
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
