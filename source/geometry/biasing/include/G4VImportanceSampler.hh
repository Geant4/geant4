#ifndef G4VImportanceSampler_hh
#define G4VImportanceSampler_hh G4VImportanceSampler_hh
#include "globals.hh"

#include "G4Nsplit_Weight.hh"


class G4VImportanceSampler{
public:
  virtual ~G4VImportanceSampler(){}
  virtual G4Nsplit_Weight Sample(G4double w) const = 0; 
};

#endif
