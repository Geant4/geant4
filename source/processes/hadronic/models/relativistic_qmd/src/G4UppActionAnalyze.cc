
#include "G4UppActionAnalyze.hh"
#include "G4VUppAnalyzer.hh"


G4UppActionAnalyze::G4UppActionAnalyze(const G4double time,
				       const G4VUppAnalyzer* aPtr, 
				       const G4UppTrackVector* sPtr)
{
  AnalyzerPtr = aPtr;
  allTracksPtr = sPtr;
  setActionTime(time);
}


G4int G4UppActionAnalyze::Perform(const G4UppTrackVector& t) const
{
  G4cout << "Starting Analyzer" << G4endl; 
  AnalyzerPtr->Analyze(*allTracksPtr);
  return 0;
}


void G4UppActionAnalyze::dump() const
{
  G4cout << "Action: Call of " << AnalyzerPtr->getName();
  G4cout << " at " << getActionTime()/fermi << G4endl;
}


