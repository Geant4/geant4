#ifndef G4VIsotopeProduction_h
#define G4VIsotopeProduction_h 1

#include "G4IsoResult.hh"
#include "G4Track.hh"
#include "G4Nucleus.hh"

class G4VIsotopeProduction
{
public:

  virtual G4IsoResult * GetIsotope(const G4Track & aTrack, const G4Nucleus & aNucleus) = 0;
  
  G4bool operator == (const G4VIsotopeProduction & aProd)
  {
    G4bool result = false;
    if(&aProd==this) result = true;
    return result;
  }
};

#endif
