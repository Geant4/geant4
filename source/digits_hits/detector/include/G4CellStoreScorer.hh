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
// $Id: G4CellStoreScorer.hh,v 1.2 2006-06-29 18:05:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4CellStoreScorer
//
// Class description:
// 
// This class invokes a G4VCellScorer according to the current
// step:
//  - pre and post step are in  a volume -> chose post gCell (cell)
//      and call ScoreInVolume
//  - cross boundary between gCells (cells)
//      call ScoreAnExtingStep on the pre cell
//      and ScoreAnEnteringStep on the post cell.
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
//

#ifndef G4CellStoreScorer_hh
#define G4CellStoreScorer_hh G4CellStoreScorer_hh


#include "G4VScorer.hh"

class G4VCellScorerStore;
class G4Step;
class G4GeometryCellStep;
class G4GeometryCell;

class G4CellStoreScorer : public G4VScorer
{
public: // with description
  explicit G4CellStoreScorer(G4VCellScorerStore &csc);
    // take a  G4VCellScorerStore created by the user to score
    // the cells contained in it. 

  virtual ~G4CellStoreScorer();
  virtual void Score(const G4Step &aStep, const G4GeometryCellStep &aPStep);
    // called to score a track for every step
 
private:
  G4VCellScorerStore &fCellScorerStore;
};



#endif

