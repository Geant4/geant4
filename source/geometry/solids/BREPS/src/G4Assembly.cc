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














































































































































