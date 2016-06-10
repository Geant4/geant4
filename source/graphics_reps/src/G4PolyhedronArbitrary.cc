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
// MODULE:                G4PolyhedronArbitrary.cc
//
// Date:                15/06/2005
// Author:                P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:                UK Ministry of Defence : RAO CRP TD Electronic Systems
// Contract:                C/MAT/N03517
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
void G4PolyhedronArbitrary::AddVertex (const G4ThreeVector& v)
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
