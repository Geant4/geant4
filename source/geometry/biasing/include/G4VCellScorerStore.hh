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
// $Id: G4VCellScorerStore.hh,v 1.7 2006-06-29 18:16:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VCellScorerStore
//
// Class description:
//
// This class describes the interface for a store for 
// G4VCellScorer objects. It is used by the class 
// G4CellStoreScorer to score the G4CellScorer objects per cell.
// 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4VCellScorerStore_hh
#define G4VCellScorerStore_hh G4VCellScorerStore_hh

class G4VCellScorer;
class G4GeometryCell;

class G4VCellScorerStore {
public: // with description 

  G4VCellScorerStore();

  virtual ~G4VCellScorerStore();

  virtual G4VCellScorer *GetCellScore(const G4GeometryCell &gCell) = 0;
    // returns a cell scorer related to the given cell 
};


#endif

