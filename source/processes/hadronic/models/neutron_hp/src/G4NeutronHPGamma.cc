// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPGamma.hh"

  G4NeutronHPGamma::G4NeutronHPGamma() 
  {
    next = NULL;
  }
  G4NeutronHPGamma::~G4NeutronHPGamma() {}

G4bool G4NeutronHPGamma::Init(ifstream & aDataFile)
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
