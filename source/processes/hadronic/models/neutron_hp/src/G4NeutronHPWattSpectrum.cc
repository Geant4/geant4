// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
#include "G4NeutronHPWattSpectrum.hh"

  G4double G4NeutronHPWattSpectrum::Sample(G4double anEnergy) 
  {
    G4double a = theApar.GetY(anEnergy)*eV;
    G4double b = theBpar.GetY(anEnergy)/eV;
    G4double result;
    G4double random, cut, max;
    max = sinh(sqrt(b*15.*a));
    do
    {
      random = G4UniformRand();
      result = -a*log(random);
      cut = G4UniformRand();
    }
    while(cut>sinh(sqrt(b*result))/max);
    return result;
  }
