// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPGamma.hh"
int G4NeutronHPGamma::instancecount = 0;
G4bool G4NeutronHPGamma::Init(G4std::ifstream & aDataFile)
{
  G4bool theResult = true;
  if(aDataFile >> levelEnergy)
  {
    aDataFile >> gammaEnergy >> probability;
    levelEnergy *= keV;
    gammaEnergy *= keV;
  }
  else
  {
    theResult=false;
  }
  return theResult;
}
