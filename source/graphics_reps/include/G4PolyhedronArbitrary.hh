// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		G4PolyhedronArbitrary.hh
//
// Date:		15/06/2005
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		UK Ministry of Defence : RAO CRP TD Electronic Systems
// Contract:		C/MAT/N03517
//
// This software is the intelectual property of QinetiQ Ltd, subject
// DEFCON 705 IPR conditions.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 13 January 2006, J Allison.  Removed unnecessary operator= functions.
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
#ifndef G4PolyhedronArbitrary_h
#define G4PolyhedronArbitrary_h 1

#include "G4Polyhedron.hh"
#include "G4ThreeVector.hh"

///////////////////////////////////////////////////////////////////////////////
//
class G4PolyhedronArbitrary : public G4Polyhedron
{
  public:
    G4PolyhedronArbitrary (const G4int nVertices, const G4int nFacets);
    virtual ~G4PolyhedronArbitrary ();
    
    void AddVertex (const G4ThreeVector v);
    void AddFacet (const G4int iv1, const G4int iv2, const G4int iv3,
      const G4int iv4=0);
    
  protected:
    G4int nVertexCount;
    G4int nFacetCount;

};
#endif
///////////////////////////////////////////////////////////////////////////////
//

