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
// $Id: G4CellStoreScorer.hh,v 1.1 2002-08-29 15:32:37 dressel Exp $
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


#include "G4VPScorer.hh"

class G4VCellScorerStore;
class G4Step;
class G4PStep;
class G4GeometryCell;

class G4CellStoreScorer : public G4VPScorer
{
public:
  G4CellStoreScorer(G4VCellScorerStore &csc);
    // take a  G4VCellScorerStore created by the user to score
    // the cells contained in it. 
  ~G4CellStoreScorer();
  void Score(const G4Step &aStep, const G4PStep &aPStep);
    // called to score a track for every step
 
private:
  G4VCellScorerStore &fCellScorerStore;
};



#endif

