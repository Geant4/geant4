#ifndef G4VImpCombiTally_hh
#define G4VImpCombiTally_hh G4VImpCombiTally_hh

#include "globals.hh"
class G4Sigma;

class G4VImpCombiTally {
public:
  G4VImpCombiTally(){}
  virtual ~G4VImpCombiTally(){}
  virtual void Tally(const G4String &rawtallyname, 
		     G4Sigma &,
		     G4double importance) = 0;
  virtual void Reset() = 0;
  virtual G4String GetName() = 0;
  virtual G4double GetValue() = 0;
  virtual G4bool HasTallied() = 0;

};

#endif
