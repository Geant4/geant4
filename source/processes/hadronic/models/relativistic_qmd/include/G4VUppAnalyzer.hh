
#ifndef G4VUPPANALYZER_H
#define G4VUPPANALYZER_H


#include "G4UppTrackVector.hh"
#include "G4UppInteraction.hh"


class G4VUppAnalyzer
{
public:

  virtual void analyze(const G4UppTrackVector& allTracks) const = 0;
  virtual void analyze(const G4UppTrackVector& allTracks, 
		       const G4UppTrackChange& aTrackChange) const = 0;

  virtual string getName() const 
    { return "Unknown Analyzer"; }

};


#endif // G4VUPPANALYZER_H
