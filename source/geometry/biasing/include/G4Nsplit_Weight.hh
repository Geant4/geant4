#ifndef G4Nsplit_Weight_hh
#define G4Nsplit_Weight_hh G4Nsplit_Weight_hh

#include "globals.hh"

struct G4Nsplit_Weight {
  G4Nsplit_Weight(G4int an, G4double aw): fN(an), fW(aw){}
  G4int fN;
  G4double fW;
};

ostream& operator<<(ostream &out, const G4Nsplit_Weight  &nw);


#endif
