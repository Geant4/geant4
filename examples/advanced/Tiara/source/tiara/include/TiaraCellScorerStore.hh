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
// $Id: TiaraCellScorerStore.hh,v 1.3 2006/06/29 15:43:22 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// Class TiaraCellScorerStore
//
// Class description:
//
// This class is a cell customized cell scorer store.
// Nevertheless the standard G4CellScorer may be created with this class
// simular to the G4CellScorerStore. In addition this class
// stores TiaraCellScorer pointers.
// This class delivers a map with the G4CellScorer objects
// to be printed with the G4ScorerTable class.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef TiaraCellScorerStore_hh
#define TiaraCellScorerStore_hh TiaraCellScorerStore_hh

#include <map>
#include "G4GeometryCell.hh"
#include "G4GeometryCellComp.hh"

#include "globals.hh"
#include "G4VCellScorerStore.hh"
#include "G4CellScorerStore.hh"

class G4CellScorer;
class IHistogram1D;
class TiaraCellScorer;

//typedef std::map<G4GeometryCell, G4CellScorer *, G4GeometryCellComp> G4MapGeometryCellCellScorer;

typedef std::map<G4GeometryCell,TiaraCellScorer  *, G4GeometryCellComp> TiaraMapGeometryCellTiaraCellScorer;

class TiaraCellScorerStore : public G4VCellScorerStore {
public:
  TiaraCellScorerStore();
  virtual ~TiaraCellScorerStore(){}
  
  virtual G4VCellScorer *GetCellScore(const G4GeometryCell &gCell);

  G4CellScorer *AddG4CellScorer(const G4GeometryCell &gCell);
  const G4MapGeometryCellCellScorer &GetMapGeometryCellCellScorer()  ;

  void AddTiaraCellScorer(TiaraCellScorer *tiaraCellScorer,
			const G4GeometryCell &gCell);
  const TiaraMapGeometryCellTiaraCellScorer &GetMapGeometryCellTiaraCellScorer() const ;
  void EndOfEventAction();

private:
  G4MapGeometryCellCellScorer fMapGeometryCellCellScorer;
  TiaraMapGeometryCellTiaraCellScorer fMapGeometryCellTiaraCellScorer;

};

#endif
