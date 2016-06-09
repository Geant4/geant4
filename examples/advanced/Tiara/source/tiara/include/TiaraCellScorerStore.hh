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
// $Id: TiaraCellScorerStore.hh,v 1.2 2003/06/18 16:40:23 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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
