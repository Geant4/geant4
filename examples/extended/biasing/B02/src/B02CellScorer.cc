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
// $Id: B02CellScorer.cc,v 1.1 2002-11-08 14:52:17 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// B02CellScorer.cc
//
// ----------------------------------------------------------------------

#include "B02CellScorer.hh"
#include "G4GeometryCell.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"


B02CellScorer::B02CellScorer(AIDA::IHistogram1D *h) :
  fHisto(h)
{}

B02CellScorer::~B02CellScorer()
{}

void B02CellScorer::ScoreAnExitingStep(const G4Step &aStep,
				       const G4GeometryCell &pre_gCell){
  fG4CellScorer.ScoreAnExitingStep(aStep, pre_gCell);
  FillHisto(aStep);
}

void B02CellScorer::ScoreAnEnteringStep(const G4Step &aStep,
					const G4GeometryCell &post_gCell){
  fG4CellScorer.ScoreAnEnteringStep(aStep, post_gCell);
  return;
}

void B02CellScorer::ScoreAnInVolumeStep(const G4Step &aStep,
					const G4GeometryCell &post_gCell){
  fG4CellScorer.ScoreAnInVolumeStep(aStep, post_gCell);
  FillHisto(aStep);
}

void B02CellScorer::FillHisto(const G4Step &aStep){

  G4StepPoint *p = 0;
  p = aStep.GetPreStepPoint();
  G4double sl = aStep.GetStepLength();
  G4double slw = sl * p->GetWeight();
  
  fHisto->fill(p->GetKineticEnergy(), slw);
}





