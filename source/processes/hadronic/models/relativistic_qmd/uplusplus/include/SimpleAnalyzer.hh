
#ifndef SIMPLEANALYZER_H
#define SIMPLEANALYZER_H

#include "G4VUppAnalyzer.hh"


class SimpleAnalyzer : public G4VUppAnalyzer
{
public:

  SimpleAnalyzer(const string& outputFile) 
    : mOutputFile(outputFile) { }
  void analyze(const G4UppTrackVector& all) const;
  void analyze(const G4UppTrackVector& all, const G4UppTrackChange& i) const
    {}
  string getName() const { return "SimpleAnalyzer"; }

private:

  const string mOutputFile;
  static int mOutputFileCnt;
};

#endif
