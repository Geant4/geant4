#ifndef G4ImpStandardTally_hh
#define G4ImpStandardTally_hh G4ImpStandardTally_hh

#include "G4VImpCombiTally.hh"
class G4StandardTally;

class G4ImpStandardTally : public G4VImpCombiTally {
public:
  G4ImpStandardTally(const G4String &tallyname, 
		     const G4String &rawtallyname,
		     const G4String &SigmaSpec);
  void Tally(const G4String &rawtallyname, 
	     G4Sigma &,
	     G4double importance);
  
  void Reset();
  G4String GetName();
  G4double GetValue();
  G4bool HasTallied();
private:
  
  G4StandardTally *fStandardTally;
  G4double fValue;


};

#endif
