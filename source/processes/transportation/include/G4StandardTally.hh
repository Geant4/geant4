#ifndef G4StandardTally_hh
#define G4StandardTally_hh G4StandardTally_hh

#include "G4VCombiTally.hh"

class G4StandardTally : public G4VCombiTally {
public:
  G4StandardTally(const G4String &tallyname, 
		  const G4String &rawtallyname,
		  const G4String &SigmaSpec);
  void Tally(const G4String &rawtallyname, 
	     G4Sigma &);

  void Reset();
  G4String GetName();
  G4double GetValue();
  G4bool HasTallied();
private:
  void Error(const G4String &m) {
    G4Exception("Error: G4StandardTally;" + fRawTallyName + ":" + m);
  }
  
  G4double ReadSigma(G4Sigma &, const G4String &sigspec);

  G4String fName;
  G4String fRawTallyName;
  G4double fValue;
  G4bool fHasTallied;
  G4String fSigmaSpec;
  
};

#endif
