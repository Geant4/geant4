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
// $Id: Tst33ScorerBuilder.cc,v 1.7 2006-06-29 22:01:20 gunter Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33ScorerBuilder.cc
//
// ----------------------------------------------------------------------

#include "Tst33ScorerBuilder.hh"
#include "G4CellScorerStore.hh"
#include "G4CellScorer.hh"
#include "Tst33VGeometry.hh"

Tst33ScorerBuilder::Tst33ScorerBuilder()
{}

Tst33ScorerBuilder::~Tst33ScorerBuilder()
{}


G4CellScorerStore *Tst33ScorerBuilder::
CreateScorer(Tst33VGeometry *samplegeo, 
	     const G4CellScorer **specialCellScorer){
  G4GeometryCell gWorldCell(samplegeo->GetWorldVolume(), -1);
  
  G4CellScorerStore *cs_store = new G4CellScorerStore();
  if (!cs_store) {
    G4Exception("Tst33ScorerBuilder::CreateScorer: new failed to create G4CellScorerStore!");
  }
  cs_store->AddCellScorer(gWorldCell);
  
  
  G4int i = 1;
  for (i=1; i <= 19; i++) {
    G4GeometryCell gCell(samplegeo->GetGeometryCell(i, ""));
    const G4CellScorer *s = cs_store->AddCellScorer(gCell);
    if (i==18) {
      *specialCellScorer = s;
    }
    if (i!=19) {
      G4GeometryCell gCellMinus(samplegeo->GetGeometryCell(i, "I1-"));
      cs_store->AddCellScorer(gCellMinus);
      G4GeometryCell gCellPlus(samplegeo->GetGeometryCell(i, "I1+"));
      cs_store->AddCellScorer(gCellPlus);
    }
  }
  return cs_store;
}
