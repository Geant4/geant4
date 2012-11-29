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
// $Id$
//
// ----------------------------------------------------------------------
// Class G4VCellScorer
//
// Class description:
//
// This is an interface for an object which does scoring for one cell
// of a geometry.
// The cell may be a physical volume or replica in the mass or a 
// parallel geometry.
//  
// Three different situation of the step with respect to the cell
// are separated: entering end exiting the cell and stepping inside
// the cell.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4VCellScorer_hh
#define G4VCellScorer_hh G4VCellScorer_hh

class G4Step;
class G4GeometryCell;

class G4VCellScorer {
public:  // with description 

  G4VCellScorer();

  virtual  ~G4VCellScorer();

  virtual void ScoreAnExitingStep(const G4Step &aStep, 
				  const G4GeometryCell &gCell) = 0;
    // to update scores related with the exiting of a cell 

  virtual void ScoreAnEnteringStep(const G4Step &aStep, 
				   const G4GeometryCell &gCell) = 0;
    //  to update scores related with the entering of a cell

  virtual void ScoreAnInVolumeStep(const G4Step &aStep, 
				   const G4GeometryCell &gCell) = 0;
    // to update scores related with the stepping inside the cell

};

#endif
