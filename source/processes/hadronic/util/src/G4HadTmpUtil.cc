#include "G4HadTmpUtil.hh"

G4int G4lrint(double ad)
  {
    return (ad>0) ? static_cast<int>(ad+.5) : static_cast<int>(ad-.5);
  }

G4int G4lint(double ad)
  {
    return (ad>0) ? static_cast<int>(ad) : static_cast<int>(ad-1.);
  }

G4int G4rint(double ad)
  {
    return (ad>0) ? static_cast<int>(ad+1) : static_cast<int>(ad);
  }
