
#ifndef G4UPPACTIONANALYZE_H
#define G4UPPACTIONANALYZE_H


#include "G4VUppAction.hh"
#include "G4VUppAnalyzer.hh"
#include "G4UppInteraction.hh"
#include "G4UppTrackVector.hh"


class G4UppActionAnalyze : public G4VUppAction
{
public:

  G4UppActionAnalyze(const G4double time,
		     const G4VUppAnalyzer* aPtr,
		     const G4UppTrackVector* sPtr);
  G4bool isValid() const { return true; }
  G4int Perform(const G4UppTrackVector& t) const;
  G4int Perform(const G4UppTrackVector& t, G4UppInteraction& i) const 
     { return 0; }
  void dump() const;
  void dump(const G4UppTrackVector& t) const { dump(); }

private:

  const G4VUppAnalyzer* AnalyzerPtr;
  const G4UppTrackVector* allTracksPtr;

};


#endif // G4UPPACTIONANALYZE_H
