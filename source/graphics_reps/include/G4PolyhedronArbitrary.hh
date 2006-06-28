//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration and of QinetiQ Ltd,  subject DEFCON 705 IPR *
// * conditions.                                                      *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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

