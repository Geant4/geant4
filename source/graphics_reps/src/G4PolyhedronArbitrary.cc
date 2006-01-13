// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		G4PolyhedronArbitrary.cc
//
// Date:		15/06/2005
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		UK Ministry of Defence : RAO CRP TD Electronic Systems
// Contract:		C/MAT/N03517
//
// This software is the intellectual property of QinetiQ Ltd, subject
// DEFCON 705 IPR conditions.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DISCLAIMER
// ----------
//
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
//
//
///////////////////////////////////////////////////////////////////////////////
//
//
#include "G4PolyhedronArbitrary.hh"
///////////////////////////////////////////////////////////////////////////////
//
G4PolyhedronArbitrary::G4PolyhedronArbitrary (const G4int nVertices,
  const G4int nFacets)
{
  AllocateMemory(nVertices, nFacets);
  nVertexCount = 0;
  nFacetCount  = 0;
}
///////////////////////////////////////////////////////////////////////////////
//
G4PolyhedronArbitrary::~G4PolyhedronArbitrary ()
{;}
///////////////////////////////////////////////////////////////////////////////
//
void G4PolyhedronArbitrary::AddVertex (const G4ThreeVector v)
{
  if (nVertexCount == nvert + 1)
  {
    G4cerr <<G4endl;
    G4cerr <<"ERROR IN G4PolyhedronArbitrary::AddVertex" <<G4endl;
    G4cerr <<"ATTEMPT TO EXCEED MAXIMUM NUMBER OF VERTICES : " << nVertexCount
           <<G4endl;
    G4cerr <<G4endl;
  }
  else
  {
    nVertexCount++;
    pV[nVertexCount] = v;
  }
}
///////////////////////////////////////////////////////////////////////////////
//
void G4PolyhedronArbitrary::AddFacet (const G4int iv1, const G4int iv2,
  const G4int iv3, const G4int iv4)
{
  if (nFacetCount == nface)
  {
    G4cerr <<G4endl;
    G4cerr <<"ERROR IN G4PolyhedronArbitrary::AddFacet" <<G4endl;
    G4cerr <<"ATTEMPT TO EXCEED MAXIMUM NUMBER OF FACETS : " << nFacetCount
           <<G4endl;
    G4cerr <<G4endl;
  }
  else if (iv1 < 1 || iv1 > nvert || iv2 < 1 || iv2 > nvert ||
            iv3 < 1 || iv3 > nvert || iv4 > nvert)
  {
    G4cerr <<G4endl;
    G4cerr <<"ERROR IN G4PolyhedronArbitrary::AddFacet" <<G4endl;
    G4cerr <<"ATTEMPT TO INDEX VERTEX NUMBER WHICH IS OUT-OF-RANGE : " <<G4endl;
    G4cerr <<G4endl;
  }
  else if (iv1 > nVertexCount || iv2 > nVertexCount || iv3 > nVertexCount ||
    iv4 > nVertexCount)
  {
    G4cerr <<G4endl;
    G4cerr <<"ERROR IN G4PolyhedronArbitrary::AddFacet" <<G4endl;
    G4cerr <<"VERTEX NEEDS TO BE DEFINED FIRST : " <<G4endl;
    G4cerr <<G4endl;
  }
  else
  {
    nFacetCount++;
    pF[nFacetCount] = G4Facet(iv1, 0, iv2, 0, iv3, 0, iv4, 0);
  }
}
///////////////////////////////////////////////////////////////////////////////
//
