// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Assembly.cc,v 1.2 2000-08-28 08:57:54 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Assembly.cc
//
// ----------------------------------------------------------------------

#include "G4Assembly.hh"

G4Assembly::G4Assembly()
{
  //  ReadSTEPFile();
  //  CopySTEPData();  
}

G4Assembly::~G4Assembly()
{
  for(G4int a=0;a<numberOfSolids;a++)
    delete placedVec[a];
}

void G4Assembly::SetPlacedVector(G4PlacedVector& pVec)
{
  numberOfSolids = pVec.entries();
  
  for(G4int a=0;a<numberOfSolids;a++)
    placedVec.append( pVec[a]);
  
}
