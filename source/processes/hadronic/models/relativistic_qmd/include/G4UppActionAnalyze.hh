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
