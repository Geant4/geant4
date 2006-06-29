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
// $Id: G4CellScorer.hh,v 1.2 2006-06-29 18:08:19 gunter Exp $
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
