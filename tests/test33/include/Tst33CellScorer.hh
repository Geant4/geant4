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
// $Id: Tst33CellScorer.hh,v 1.1 2007-06-22 12:47:15 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class Tst33CellScorer
//
// Class description:
//
// This class is an example how to build a customized cell scorer
// derived from G4VCellScorer that also uses the G4CellScorer.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33CellScorer_hh
#define Tst33CellScorer_hh Tst33CellScorer_hh

#include "G4VCellScorer.hh"
#include "G4CellScoreComposer.hh"
#include "G4CellScorer.hh"

#include <vector>

class Tst33CellScorer : public G4VCellScorer {
public:  
  Tst33CellScorer();
  virtual ~Tst33CellScorer();
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

  G4CellScorer fG4CellScorer;
};



#endif
