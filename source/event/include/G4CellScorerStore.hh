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
// $Id: G4CellScorerStore.hh,v 1.1 2002-10-28 10:06:00 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4CellScorerStore
//
// Class description:
//
// This class stores pointers to G4CellScorer. Objects of G4CellScorer for
// each cell may be created the first time a cell is hit: therefore 
// SetAutoScorerCreate() has to be called before the run starts.
// Alternative the store may be filled with pointers to  G4CellScorer
// before the run. In this case only the cells for which a G4CellScorer
// exists are scored.
// 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4CellScorerStore_hh
#define G4CellScorerStore_hh G4CellScorerStore_hh

#include "g4std/map"
#include "G4GeometryCell.hh"

#include "globals.hh"
#include "G4VCellScorerStore.hh"
#include "G4GeometryCellComp.hh"

class G4CellScorer;

typedef G4std::map<G4GeometryCell, G4CellScorer *, G4GeometryCellComp> G4MapGeometryCellCellScorer;

class G4CellScorerStore : public G4VCellScorerStore {
public: // with description

  G4CellScorerStore();

  virtual ~G4CellScorerStore();
  
  virtual G4VCellScorer *GetCellScore(const G4GeometryCell &gCell);
    // return the current cell scorer
 
  void SetAutoScorerCreate();
    // if this function is called cell scorers will be created 
    // automatically for a cell the first time a track hits the cell 

  G4CellScorer *AddCellScorer(const G4GeometryCell &gCell);
    // add a G4CellScorer for the given cell, needed if 
    // SetAutoScorerCreate is not called. This way the user
    // may specify to score only for the given cells.


  G4CellScorer *AddCellScorer(G4VPhysicalVolume &vol, G4int repnum = 0);
    // same functionallity as above only that the cell is specified
    // by the physical volume and a replica number


  const G4MapGeometryCellCellScorer &GetMapGeometryCellCellScorer() const ;
    // return the map of GeometryCell to the G4CellScorer

  void DeleteAllScorers();

private:
  G4MapGeometryCellCellScorer fMapGeometryCellCellScorer;
  G4bool fAutoCreate;
};

#endif
