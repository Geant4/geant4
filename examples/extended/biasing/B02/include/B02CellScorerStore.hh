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
// $Id: B02CellScorerStore.hh,v 1.1 2002-11-08 14:52:16 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class B02CellScorerStore
//
// Class description:
//
// This class is a cell customized cell scorer store.
// Nevertheless the standard G4CellScorer may be created with this class
// simular to the G4CellScorerStore. In addition this class
// stores B02CellScorer pointers.
// This class delivers a map with the G4CellScorer objects
// to be printed with the G4ScorerTable class.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef B02CellScorerStore_hh
#define B02CellScorerStore_hh B02CellScorerStore_hh

#include "g4std/map"
#include "G4GeometryCell.hh"
#include "G4GeometryCellComp.hh"

#include "globals.hh"
#include "G4VCellScorerStore.hh"

class G4CellScorer;
class B02CellScorer;
class IHistogram1D;

typedef G4std::map<G4GeometryCell, G4CellScorer *, G4GeometryCellComp> G4MapGeometryCellCellScorer;

typedef G4std::map<G4GeometryCell,B02CellScorer  *, G4GeometryCellComp> B02MapGeometryCellB02CellScorer;

class B02CellScorerStore : public G4VCellScorerStore {
public:
  B02CellScorerStore();
  virtual ~B02CellScorerStore(){}
  
  virtual G4VCellScorer *GetCellScore(const G4GeometryCell &gCell);

  G4CellScorer *AddG4CellScorer(const G4GeometryCell &gCell);
  const G4MapGeometryCellCellScorer &GetMapGeometryCellCellScorer()  ;

  void AddB02CellScorer(B02CellScorer *b02cellScorer,
			const G4GeometryCell &gCell);
  const B02MapGeometryCellB02CellScorer &GetMapGeometryCellB02CellScorer() const ;
  
private:
  G4MapGeometryCellCellScorer fMapGeometryCellCellScorer;
  B02MapGeometryCellB02CellScorer fMapGeometryCellB02CellScorer;

};

#endif
