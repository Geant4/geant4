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
// $Id: Tst33CellScorerStore.hh,v 1.1 2007-06-22 12:47:15 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class Tst33CellScorerStore
//
// Class description:
//
// This class is a cell customized cell scorer store.
// Nevertheless the standard G4CellScorer may be created with this class
// simular to the G4CellScorerStore. In addition this class
// stores Tst33CellScorer pointers.
// This class delivers a map with the G4CellScorer objects
// to be printed with the G4ScorerTable class.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33CellScorerStore_hh
#define Tst33CellScorerStore_hh Tst33CellScorerStore_hh

#include <map>
#include "G4GeometryCell.hh"
#include "G4GeometryCellComp.hh"

#include "globals.hh"
#include "G4VCellScorerStore.hh"

class G4CellScorer;
class Tst33CellScorer;
class IHistogram1D;

typedef std::map<G4GeometryCell, G4CellScorer *, G4GeometryCellComp> G4MapGeometryCellCellScorer;

typedef std::map<G4GeometryCell,Tst33CellScorer  *, G4GeometryCellComp> Tst33MapGeometryCellTst33CellScorer;

class Tst33CellScorerStore : public G4VCellScorerStore {
public:
  Tst33CellScorerStore();
  virtual ~Tst33CellScorerStore(){}
  
  virtual G4VCellScorer *GetCellScore(const G4GeometryCell &gCell);

  G4CellScorer *AddG4CellScorer(const G4GeometryCell &gCell);
  const G4MapGeometryCellCellScorer &GetMapGeometryCellCellScorer()  ;

  void AddTst33CellScorer(Tst33CellScorer *b02cellScorer,
			const G4GeometryCell &gCell);
  const Tst33MapGeometryCellTst33CellScorer &GetMapGeometryCellTst33CellScorer() const ;
  
private:
  G4MapGeometryCellCellScorer fMapGeometryCellCellScorer;
  Tst33MapGeometryCellTst33CellScorer fMapGeometryCellTst33CellScorer;

};

#endif
