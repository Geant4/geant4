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
// $Id: B02CellScorer.hh,v 1.1 2002-11-08 14:52:16 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class B02CellScorer
//
// Class description:
//
// This class is an example how to build a customized cell scorer
// derived from G4VCellScorer that also uses the G4CellScorer.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef B02CellScorer_hh
#define B02CellScorer_hh B02CellScorer_hh

#include "G4VCellScorer.hh"
#include "G4CellScoreComposer.hh"
#include "G4CellScorer.hh"

#include "g4std/vector"

#include "AIDA/IHistogram1D.h"

class B02CellScorer : public G4VCellScorer {
public:  
  B02CellScorer(AIDA::IHistogram1D *h);
  virtual ~B02CellScorer();
  virtual void ScoreAnExitingStep(const G4Step &aStep, 
				  const G4GeometryCell &gCell);
  virtual void ScoreAnEnteringStep(const G4Step &aStep, 
				   const G4GeometryCell &gCell);
  virtual void ScoreAnInVolumeStep(const G4Step &aStep, 
				   const G4GeometryCell &gCell);

  G4CellScorer &GetG4CellScorer(){
    return fG4CellScorer;
  }
private:
  void FillHisto(const G4Step &aStep);

  G4CellScorer fG4CellScorer;
  AIDA::IHistogram1D *fHisto;
};



#endif
