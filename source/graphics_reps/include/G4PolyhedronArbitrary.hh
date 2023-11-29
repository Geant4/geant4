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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:                G4PolyhedronArbitrary.hh
//
// Date:                15/06/2005
// Author:                P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:                UK Ministry of Defence : RAO CRP TD Electronic Systems
// Contract:                C/MAT/N03517
//
// This software is the intelectual property of QinetiQ Ltd, subject
// DEFCON 705 IPR conditions.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 11th November 2011, J Allison.  Added private copy constructor and
//    assignment operator added to satisfy Coverity.
// 27th July 2011, J Allison.  Added SetReferences and InvertFacets.
//   SetReferences is necessary at the end to complete the polyhedron.
//   It particularly matters if the polyhedron suffers subsequent
//   boolean operations.
//   InvertFacets can be useful.
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
    ~G4PolyhedronArbitrary () override;
    // Private copy constructor and assignment operator added to satisfy
    // Coverity - JA 11/11/11.
    G4PolyhedronArbitrary(const G4PolyhedronArbitrary&) = delete;
    G4PolyhedronArbitrary& operator= (const G4PolyhedronArbitrary&) = delete;    

    void AddVertex (const G4ThreeVector& v);
    void AddFacet (const G4int iv1, const G4int iv2, const G4int iv3,
      const G4int iv4=0);
    
    // Call this after all vertices and facets have been added.
    void SetReferences() {HepPolyhedron::SetReferences();}

    // Can be useful.
    void InvertFacets() {HepPolyhedron::InvertFacets();}

  protected:
    G4int nVertexCount;
    G4int nFacetCount;

};
#endif
///////////////////////////////////////////////////////////////////////////////
//

