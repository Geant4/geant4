#ifndef G4ImportanceTally_hh
#define G4ImportanceTally_hh

#include "G4VImpCombiTally.hh"

class G4ImportanceTally : public G4VImpCombiTally{
public:
  G4ImportanceTally(const G4String &tallyname);
  void Tally(const G4String &rawtallyname, 
	     G4Sigma &,
	     G4double importance);
  void Reset();
  G4String GetName();
  G4double GetValue();
  G4bool HasTallied();
private:
  G4String fName;
  G4double fImportance;
  G4bool fHasTallied;
  
};
#endif
