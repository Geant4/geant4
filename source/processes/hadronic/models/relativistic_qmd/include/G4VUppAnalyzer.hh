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
