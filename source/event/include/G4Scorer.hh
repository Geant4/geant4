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
//
// $Id: G4Scorer.hh,v 1.1 2002-10-28 10:06:01 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Scorer
//
// Class description:
//
// This is a simple scorer using the G4CellScorerStore and 
// G4CellStoreScorer.
// 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4Scorer_hh
#define G4Scorer_hh G4Scorer_hh

#include "G4VScorer.hh"
#include "G4CellScorerStore.hh"
#include "G4CellStoreScorer.hh"

class G4Scorer : public G4VScorer {
public: // with description
  G4Scorer();

  virtual ~G4Scorer();

  virtual void Score(const G4Step &aStep, const G4GeometryCellStep &aPStep);
    // called to score a track for every step

  const G4MapGeometryCellCellScorer &GetMapGeometryCellCellScorer() const ;
    // returns a reference to the store containing the scores

private:
  G4Scorer(const G4Scorer &);
  G4Scorer &operator=(const G4Scorer &);
  G4CellScorerStore fCellScorerStore;  
  G4CellStoreScorer fCellStoreScorer;
};

#endif
