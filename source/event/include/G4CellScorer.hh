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
// $Id: G4CellScorer.hh,v 1.1 2002-10-28 10:06:00 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4CellScorer
//
// Class description:
//
// This class does standard scoring for a cell. The cell may be a 
// physical volume or replica in the mass or a parallel geometry.
// One object of this class per cell is used.
// See intercoms/include/G4VCellScorer.hh
// 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4CellScorer_hh
#define G4CellScorer_hh G4CellScorer_hh

#include "G4VCellScorer.hh"
#include "G4TrackLogger.hh"
#include "G4CellScoreComposer.hh"
#include "G4CellScoreValues.hh"

class G4CellScorer : public G4VCellScorer {
public:  // with description

  G4CellScorer();

  virtual ~G4CellScorer();

  virtual void ScoreAnExitingStep(const G4Step &aStep, 
				  const G4GeometryCell &gCell);
   // to update scores related with the exiting of a cell 

  virtual void ScoreAnEnteringStep(const G4Step &aStep, 
				   const G4GeometryCell &gCell);
    //  to update scores related with the entering of a cell

  virtual void ScoreAnInVolumeStep(const G4Step &aStep, 
				   const G4GeometryCell &gCell);
    // to update scores related with the stepping inside the cell

  const G4CellScoreComposer &GetCellScoreComposer() const;
    // return the G4CellScoreComposer 

  const G4CellScoreValues &GetCellScoreValues() const;
    // return scores in G4CellScoreValues

private:
  void ScorePopulation(const G4GeometryCell &post_gCell, 
		       const G4Step &aStep);
  G4CellScoreComposer fCellScoreComposer;
  G4TrackLogger fTrackLogger;
};



#endif
