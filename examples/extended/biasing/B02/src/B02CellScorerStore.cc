//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file biasing/B02/src/B02CellScorerStore.cc
/// \brief Implementation of the B02CellScorerStore class
//
//
// $Id$
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
}

void B02CellScorerStore::
AddB02CellScorer(B02CellScorer *b02scorer,
                 const G4GeometryCell &g) {
  fMapGeometryCellB02CellScorer[g] = b02scorer;
}


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
