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
// $Id: Tst33CellScorerStore.cc,v 1.2 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33CellScorerStore.cc
//
// ----------------------------------------------------------------------

#include "Tst33CellScorerStore.hh"

#include "G4CellScorer.hh"
#include "Tst33CellScorer.hh"

Tst33CellScorerStore::Tst33CellScorerStore()
{}


G4CellScorer *Tst33CellScorerStore::
AddG4CellScorer(const G4GeometryCell &g) {
  return fMapGeometryCellCellScorer[g] = 
    new G4CellScorer;
}

void Tst33CellScorerStore::
AddTst33CellScorer(Tst33CellScorer *b02scorer,
		 const G4GeometryCell &g) {
  fMapGeometryCellTst33CellScorer[g] = b02scorer;
}


const G4MapGeometryCellCellScorer &Tst33CellScorerStore::GetMapGeometryCellCellScorer()  {      
  for (Tst33MapGeometryCellTst33CellScorer::iterator it = fMapGeometryCellTst33CellScorer.begin();
       it != fMapGeometryCellTst33CellScorer.end(); it++) {
    fMapGeometryCellCellScorer[(*it).first] = 
      &((*it).second->GetG4CellScorer());
  }
  return fMapGeometryCellCellScorer;
}

const Tst33MapGeometryCellTst33CellScorer &Tst33CellScorerStore::GetMapGeometryCellTst33CellScorer() const {
  return fMapGeometryCellTst33CellScorer;
}

G4VCellScorer *Tst33CellScorerStore::
GetCellScore(const G4GeometryCell &gCell){
  Tst33MapGeometryCellTst33CellScorer::iterator itb08 = 
    fMapGeometryCellTst33CellScorer.find(gCell);
  if (itb08 != fMapGeometryCellTst33CellScorer.end()) {
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
