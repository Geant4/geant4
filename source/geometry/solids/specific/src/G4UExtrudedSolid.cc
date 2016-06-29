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
//
// $Id:$
//
// 
// Implementation of G4UExtrudedSolid wrapper class
// --------------------------------------------------------------------

#include "G4ExtrudedSolid.hh"
#include "G4UExtrudedSolid.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4PolyhedronArbitrary.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructors
//
G4UExtrudedSolid::G4UExtrudedSolid(const G4String&          name,
                                   std::vector<G4TwoVector> polygon,
                                   std::vector<ZSection>    zsections)
  : G4USolid(name, new UExtrudedSolid())
{
  GetShape()->SetName(name);
  std::vector<UVector2> pvec;
  for (unsigned int i=0; i<polygon.size(); ++i)
  {
    pvec.push_back(UVector2(polygon[i].x(), polygon[i].y()));
  }
  std::vector<UExtrudedSolid::ZSection> svec;
  for (unsigned int i=0; i<zsections.size(); ++i)
  {
    ZSection sec = zsections[i];
    svec.push_back(UExtrudedSolid::ZSection(sec.fZ,
                   UVector2(sec.fOffset.x(), sec.fOffset.y()), sec.fScale));
  }
  GetShape()->Initialise(pvec, svec);
}


G4UExtrudedSolid::G4UExtrudedSolid(const G4String&          name,
                                   std::vector<G4TwoVector> polygon,
                                   G4double                 halfZ,
                                   G4TwoVector off1, G4double scale1,
                                   G4TwoVector off2, G4double scale2)
  : G4USolid(name, new UExtrudedSolid())
{ 
  GetShape()->SetName(name);
  std::vector<UVector2> pvec;
  for (unsigned int i=0; i<polygon.size(); ++i)
  {
    pvec.push_back(UVector2(polygon[i].x(), polygon[i].y()));
  }
  GetShape()->Initialise(pvec, halfZ, UVector2(off1.x(), off1.y()), scale1,
                                      UVector2(off2.x(), off2.y()), scale2);
}

////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UExtrudedSolid::G4UExtrudedSolid(__void__& a)
  : G4USolid(a)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UExtrudedSolid::~G4UExtrudedSolid()
{
}


//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UExtrudedSolid::G4UExtrudedSolid(const G4UExtrudedSolid &source)
  : G4USolid(source)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UExtrudedSolid&
G4UExtrudedSolid::operator=(const G4UExtrudedSolid &source)
{
  if (this == &source) return *this;
  
  G4USolid::operator=( source );
  
  return *this;
}


//////////////////////////////////////////////////////////////////////////
//
// Accessors

G4int G4UExtrudedSolid::GetNofVertices() const
{
  return GetShape()->GetNofVertices();
}
G4TwoVector G4UExtrudedSolid::GetVertex(G4int i) const
{
  UVector2 v = GetShape()->GetVertex(i);
  return G4TwoVector(v.x, v.y);
}
std::vector<G4TwoVector> G4UExtrudedSolid::GetPolygon() const
{
  std::vector<UVector2> pol = GetShape()->GetPolygon();
  std::vector<G4TwoVector> v;
  for (unsigned int i=0; i<pol.size(); ++i)
  {
    v.push_back(G4TwoVector(pol[i].x, pol[i].y));
  }
  return v;
}
G4int G4UExtrudedSolid::GetNofZSections() const
{
  return GetShape()->GetNofZSections();
}
G4UExtrudedSolid::ZSection G4UExtrudedSolid::GetZSection(G4int i) const
{
  return ZSection(GetShape()->GetZSection(i));
}
std::vector<G4UExtrudedSolid::ZSection> G4UExtrudedSolid::GetZSections() const
{
  std::vector<UExtrudedSolid::ZSection> sv = GetShape()->GetZSections();
  std::vector<G4UExtrudedSolid::ZSection> vec;
  for (unsigned int i=0; i<sv.size(); ++i)
  {
    vec.push_back(ZSection(sv[i]));
  }
  return vec;
}


///////////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron()
//
G4Polyhedron* G4UExtrudedSolid::CreatePolyhedron () const
{
  G4int nFacets = GetShape()->GetNumberOfFacets();
  G4int nVertices = 0;
  for (G4int l = 0; l<nFacets; ++l)  // compute total number of vertices first
  {
    VUFacet* facet = GetShape()->GetFacet(l);
    G4int n = facet->GetNumberOfVertices();
    nVertices += n;
  }

  G4PolyhedronArbitrary *polyhedron =
    new G4PolyhedronArbitrary (nVertices,nFacets);

  for (G4int i = 0; i<nFacets; ++i)
  {
    VUFacet* facet = GetShape()->GetFacet(i);
    G4int v[4];
    G4int n = facet->GetNumberOfVertices();
    for (G4int m = 0; m<n; ++m)
    {
      UVector3 vtx = facet->GetVertex(m);
      polyhedron->AddVertex(G4ThreeVector(vtx.x(), vtx.y(), vtx.z()));
    }
    if (n > 4) n = 4;
    else if (n == 3) v[3] = 0;
    for (G4int j=0; j<n; ++j)
    {
      G4int k = facet->GetVertexIndex(j);
      v[j] = k+1;
    }
    polyhedron->AddFacet(v[0],v[1],v[2],v[3]);
  }
  polyhedron->SetReferences();  

  return (G4Polyhedron*) polyhedron;
}

#endif  // G4GEOM_USE_USOLIDS
