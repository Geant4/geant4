#include "G4tgbMaterial.hh"

#include "G4tgrMaterial.hh"


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbMaterial::G4tgbMaterial( G4tgrMaterial* hg )
{
  theTgrMate = hg;
  theG4Mate = 0;
}

