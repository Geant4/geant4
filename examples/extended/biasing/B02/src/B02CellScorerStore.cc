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
// $Id: B02CellScorerStore.cc,v 1.1 2002-11-08 14:52:17 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// B02CellScorerStore.cc
//
// ----------------------------------------------------------------------

#include "B02CellScorerStore.hh"

#include "G4CellScorer.hh"
#include "B02CellScorer.hh"

B02CellScorerStore::B02CellScorerStore()
{}


G4CellScorer *B02CellScorerStore::
AddG4CellScorer(const G4GeometryCell &g) {
  return fMapGeometryCellCellScorer[g] = 
    new G4CellScorer;
};

void B02CellScorerStore::
AddB02CellScorer(B02CellScorer *b02scorer,
		 const G4GeometryCell &g) {
  fMapGeometryCellB02CellScorer[g] = b02scorer;
};


const G4MapGeometryCellCellScorer &B02CellScorerStore::GetMapGeometryCellCellScorer()  {      
  for (B02MapGeometryCellB02CellScorer::iterator it = fMapGeometryCellB02CellScorer.begin();
       it != fMapGeometryCellB02CellScorer.end(); it++) {
    fMapGeometryCellCellScorer[(*it).first] = 
      &((*it).second->GetG4CellScorer());
  }
  return fMapGeometryCellCellScorer;
}

const B02MapGeometryCellB02CellScorer &B02CellScorerStore::GetMapGeometryCellB02CellScorer() const {
  return fMapGeometryCellB02CellScorer;
}

G4VCellScorer *B02CellScorerStore::
GetCellScore(const G4GeometryCell &gCell){
  B02MapGeometryCellB02CellScorer::iterator itb08 = 
    fMapGeometryCellB02CellScorer.find(gCell);
  if (itb08 != fMapGeometryCellB02CellScorer.end()) {
    return (*itb08).second;
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
