#ifndef G4VIsotopeProduction_h
#define G4VIsotopeProduction_h 1

#include "G4IsoResult.hh"
#include "G4Track.hh"
#include "G4Nucleus.hh"

// Class Description
// This is the class you inherit from, if you want to implement your special
// isotope production model, based on your production cross-sections.
// Registering it with the corresponding process, the particle flux from the
// transport models will be fed into the production model, to retrieve improved
// isotope production information.
// Class Description - End
class G4VIsotopeProduction
{
public:// With Description

  // This is the interface to implement for isotope production models.
  virtual G4IsoResult * GetIsotope(const G4Track & aTrack, const G4Nucleus & aNucleus) = 0;
  
public:// Without Description
  G4bool operator == (const G4VIsotopeProduction & aProd)
  {
    G4bool result = false;
    if(&aProd==this) result = true;
    return result;
  }
};

#endif
