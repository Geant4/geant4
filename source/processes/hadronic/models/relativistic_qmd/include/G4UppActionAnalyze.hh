
#ifndef G4UPPACTIONANALYZE_H
#define G4UPPACTIONANALYZE_H


#include "G4VUppAction.hh"
#include "G4VUppAnalyzer.hh"
#include "G4UppInteraction.hh"
#include "G4UppTrackVector.hh"


class G4UppActionAnalyze : public G4VUppAction
{
public:

  G4UppActionAnalyze(const G4double analyzeTime,
		     const G4VUppAnalyzer& anAnalzer);

  G4bool isValid() const 
    { return true; }

  G4UppTrackChange* perform(const G4UppTrackVector& allTracks) const;

  void dump() const;

private:

  const G4VUppAnalyzer* analyzerPtr;

};


#endif // G4UPPACTIONANALYZE_H
