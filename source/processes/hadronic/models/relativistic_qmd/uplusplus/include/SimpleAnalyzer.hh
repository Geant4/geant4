//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

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
