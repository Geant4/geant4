// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Assembly.cc,v 1.4 2001-04-20 19:55:26 gcosmo Exp $
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
  G4PlacedSolid* a = 0;
  
  // Remove placedVec and delete all its contents
  while (placedVec.size()>0)
  {
    a = placedVec.back();
    placedVec.pop_back();
    for (G4PlacedVector::iterator i=placedVec.begin(); i!=placedVec.end(); i++)
    {
      if (*i==a)
      {
	placedVec.erase(i);
	i--;
      }
    } 
    if ( a )  delete a;    
  } 
}

void G4Assembly::SetPlacedVector(G4PlacedVector& pVec)
{
  numberOfSolids = pVec.size();
  
  for(G4int a=0;a<numberOfSolids;a++)
    placedVec.push_back( pVec[a]);
  
}
