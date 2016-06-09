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
//
// $Id: TiaraCellScorerStore.cc,v 1.1.1.2.2.1 2008/09/03 08:40:16 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-01-patch-03 $
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
}

void TiaraCellScorerStore::
AddTiaraCellScorer(TiaraCellScorer *tiaraScorer,
		 const G4GeometryCell &g) {
  fMapGeometryCellTiaraCellScorer[g] = tiaraScorer;
}

void TiaraCellScorerStore::EndOfEventAction() {
  for (TiaraMapGeometryCellTiaraCellScorer::iterator it = fMapGeometryCellTiaraCellScorer.begin();
       it != fMapGeometryCellTiaraCellScorer.end(); it++) {
    (*it).second->EndOfEventAction();
  }
}

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
