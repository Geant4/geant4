
#include "G4UppActionAnalyze.hh"
#include "G4VUppAnalyzer.hh"


G4UppActionAnalyze::G4UppActionAnalyze(const G4double analyzeTime,
				       const G4VUppAnalyzer& anAnalyzer)
{
  analyzerPtr = &anAnalyzer;
  setActionTime(analyzeTime);
}


G4UppTrackChange* G4UppActionAnalyze::perform(const G4UppTrackVector& allTracks) const
{
  // G4cout << "Starting Analyzer" << G4endl; 
  analyzerPtr->analyze(allTracks);
  return NULL;
}


void G4UppActionAnalyze::dump() const
{
  G4cout << "Action: ANALYSE (" << analyzerPtr->getName();
  G4cout << ") at " << getActionTime()*c_light/fermi << " fm/c" << G4endl;
}


