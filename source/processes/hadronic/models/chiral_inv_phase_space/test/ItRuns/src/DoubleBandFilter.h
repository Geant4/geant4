#ifndef DoubleBandFilter_h
#define DoubleBandFilter_h

#include "Analysis/src/VFilter.h"
#include "globals.hh"

class DoubleBandFilter : public TVANAFilter<G4double>
{
  public:
  
  DoubleBandFilter(G4double higherThan, G4double lowerEquals, G4String aName)
   : TVANAFilter<G4double>(aName)
  {
    theLow = higherThan;
    theHigh = lowerEquals;
  }
  
  G4bool Accept(G4double & anInput)
  {
    if(anInput<=theHigh && anInput>theLow) return true;
    return false;
  }
  
  private:
  
  G4double theLow;
  G4double theHigh;
};

#endif
