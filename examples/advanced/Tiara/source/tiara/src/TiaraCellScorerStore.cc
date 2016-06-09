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
// $Id: TiaraCellScorerStore.cc,v 1.1.1.1 2003/06/12 13:08:25 dressel Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// TiaraCellScorerStore.cc
//
// ----------------------------------------------------------------------

#include "TiaraCellScorerStore.hh"

#include "G4CellScorer.hh"
#include "TiaraCellScorer.hh"


TiaraCellScorerStore::TiaraCellScorerStore()
{}


G4CellScorer *TiaraCellScorerStore::
AddG4CellScorer(const G4GeometryCell &g) {
  return fMapGeometryCellCellScorer[g] = 
    new G4CellScorer;
};

void TiaraCellScorerStore::
AddTiaraCellScorer(TiaraCellScorer *tiaraScorer,
		 const G4GeometryCell &g) {
  fMapGeometryCellTiaraCellScorer[g] = tiaraScorer;
};

void TiaraCellScorerStore::EndOfEventAction() {
  for (TiaraMapGeometryCellTiaraCellScorer::iterator it = fMapGeometryCellTiaraCellScorer.begin();
       it != fMapGeometryCellTiaraCellScorer.end(); it++) {
    (*it).second->EndOfEventAction();
  }
};

const G4MapGeometryCellCellScorer &TiaraCellScorerStore::GetMapGeometryCellCellScorer()  {      
  for (TiaraMapGeometryCellTiaraCellScorer::iterator it = fMapGeometryCellTiaraCellScorer.begin();
       it != fMapGeometryCellTiaraCellScorer.end(); it++) {
    fMapGeometryCellCellScorer[(*it).first] = 
      new G4CellScorer((*it).second->GetG4CellScorer());
  }
  return fMapGeometryCellCellScorer;
}

const TiaraMapGeometryCellTiaraCellScorer &TiaraCellScorerStore::GetMapGeometryCellTiaraCellScorer() const {
  return fMapGeometryCellTiaraCellScorer;
}

G4VCellScorer *TiaraCellScorerStore::
GetCellScore(const G4GeometryCell &gCell){
  TiaraMapGeometryCellTiaraCellScorer::iterator itTiara = 
    fMapGeometryCellTiaraCellScorer.find(gCell);
  if (itTiara != fMapGeometryCellTiaraCellScorer.end()) {
    return (*itTiara).second;
  }
  else {
    G4MapGeometryCellCellScorer::iterator it = fMapGeometryCellCellScorer.find(gCell);
    if (it != fMapGeometryCellCellScorer.end()) {
      return (*it).second;
    }
    else {
      return 0;
    }
  }
}
