#ifndef G4VAnnihilationCrossSection_h
#define G4VAnnihilationCrossSection_h

#include "globals.hh"

class G4VAnnihilationCrossSection
{
  public:
    virtual G4bool InCharge(G4int aCode, G4int bCode) = 0;
    virtual G4double GetXsec(G4double s) = 0;
};

#endif
