
#ifndef G4VUPPANALYZER_H
#define G4VUPPANALYZER_H


#include "G4UppTrackVector.hh"
#include "G4UppInteraction.hh"


class G4VUppAnalyzer
{
public:

  virtual void Analyze(const G4UppTrackVector& s) const = 0;
  virtual void Analyze(const G4UppTrackVector& s, const G4UppInteraction& i) const = 0;
  virtual string getName() const { return "Unknown Analyzer"; }

};


#endif // G4VUPPANALYZER_H
