#ifndef G4VCombiTally_hh
#define G4VCombiTally_hh G4VCombiTally_hh

#include "globals.hh"
class G4Sigma;

class G4VCombiTally {
public:
  G4VCombiTally(){}
  virtual ~G4VCombiTally(){}
  virtual void Tally(const G4String &rawtallyname, 
		     G4Sigma &) = 0;
  virtual void Reset() = 0;
  virtual G4String GetName() = 0;
  virtual G4double GetValue() = 0;
  virtual G4bool HasTallied() = 0;

};

#endif
