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
// $Id$
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4CellScorerStore.cc
//
// ----------------------------------------------------------------------

#include "G4CellScorerStore.hh"

#include "G4CellScorer.hh"

G4CellScorerStore::G4CellScorerStore() :
  fAutoCreate(false)
{
  G4cout << "--------------------------------------------------------" << G4endl
         << "WARNING: Class  <G4CellScorerStore>  is  now obsolete |" << G4endl
         << "         and will be removed starting from next Geant4 |" << G4endl
         << "         major release.  Please, consider switching to |" << G4endl
         << "         general purpose scoring functionality.        |" << G4endl
         << "--------------------------------------------------------"
         << G4endl;
}

G4CellScorerStore::~G4CellScorerStore()
{}

void G4CellScorerStore::SetAutoScorerCreate(){
  fAutoCreate = true;
}

G4CellScorer *G4CellScorerStore::
AddCellScorer(G4VPhysicalVolume &vol, G4int repnum) {
  G4CellScorer *cs = 0;
  cs = new G4CellScorer;
  if (!cs) {
    G4Exception("G4CellScorerStore::AddCellScorer","Event0801",FatalException,
    "failed to create G4CellScorer!");
  }  
  fMapGeometryCellCellScorer[G4GeometryCell(vol, repnum)] = cs;
  return cs;
}

G4CellScorer *G4CellScorerStore::
AddCellScorer(const G4GeometryCell &gCell) {
  G4CellScorer *cs = 0;
  cs = new G4CellScorer;
  if (!cs) {
    G4Exception("G4CellScorerStore::AddCellScorer","Event0801",FatalException,
    "failed to create G4CellScorer!");
  }  
  fMapGeometryCellCellScorer[gCell] = cs;
  return cs;
}

const G4MapGeometryCellCellScorer &G4CellScorerStore::GetMapGeometryCellCellScorer() const {
  return fMapGeometryCellCellScorer;
}

G4VCellScorer *G4CellScorerStore::
GetCellScore(const G4GeometryCell &gCell){
  G4VCellScorer *cs=0;
  G4MapGeometryCellCellScorer::iterator it = fMapGeometryCellCellScorer.find(gCell);
  if (it != fMapGeometryCellCellScorer.end()) {
    cs = (*it).second;
  }
  else {
    if (fAutoCreate) {
      cs =  AddCellScorer(gCell);
    }
  }
  return cs;
}

void G4CellScorerStore::DeleteAllScorers(){
  for(G4MapGeometryCellCellScorer::iterator it = 
	fMapGeometryCellCellScorer.begin(); 
      it!=fMapGeometryCellCellScorer.end();
      ++it){
    delete it->second;
  }
  fMapGeometryCellCellScorer.clear();
}
