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
// $Id: G4Scorer.cc,v 1.1 2002-10-28 10:06:01 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Scorer.cc
//
// ----------------------------------------------------------------------

#include "G4Scorer.hh"

G4Scorer::G4Scorer() 
  :
  fCellStoreScorer(fCellScorerStore)
{
  fCellScorerStore.SetAutoScorerCreate();
}

G4Scorer::~G4Scorer()
{
}

void G4Scorer::Score(const G4Step &aStep, const G4GeometryCellStep &aPStep){
  fCellStoreScorer.Score(aStep, aPStep);
}

const G4MapGeometryCellCellScorer &G4Scorer::GetMapGeometryCellCellScorer() const {
  return fCellScorerStore.GetMapGeometryCellCellScorer();
}

