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


