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
// $Id: G4StandardScorer.hh,v 1.1 2002-07-10 15:51:04 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4StandardScorer
//
// Class description:
// 
// This class steers the standard scoring.
// 
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
//

#ifndef G4StandardScorer_hh
#define G4StandardScorer_hh G4StandardScorer_hh

#include "g4std/iostream"
#include "g4std/map"

#include "G4VPScorer.hh"
#include "G4MapPtkStandardCellScorer.hh"
#include "G4TrackLogger.hh"

class G4Step;
class G4PStep;
class G4VIStore;

class G4StandardScorer : public G4VPScorer
{
public:
  G4StandardScorer();
  ~G4StandardScorer();
  void Score(const G4Step &aStep, const G4PStep &aPStep);
    // called to score a track for every step
 
  const G4MapPtkStandardCellScorer &GetMapPtkStandardCellScorer() const 
  { return fPtkScores; }
    // returns a the G4MapPtkStandardCellScorer cotaining
    // G4StandardCellScorer for every cell

  typedef G4std::map<G4PTouchableKey, G4TrackLogger ,G4PTkComp > MapPtkTrackLogger; 
    // used to avoid counting reantrance of tracks in the population

private:
  G4MapPtkStandardCellScorer fPtkScores;
  MapPtkTrackLogger fMapPtkTrackLogger;
};

G4std::ostream& operator<<(G4std::ostream &out, 
			   const G4StandardScorer &ps);


#endif

