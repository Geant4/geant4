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
// $Id: G4CellScorer.cc,v 1.2 2002-11-23 12:30:40 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4CellScorer.cc
//
// ----------------------------------------------------------------------

#include "G4CellScorer.hh"
#include "G4GeometryCell.hh"
#include "G4Track.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

G4CellScorer::G4CellScorer()
{}

G4CellScorer::~G4CellScorer()
{}

void G4CellScorer::ScoreAnExitingStep(const G4Step &aStep,
				      const G4GeometryCell &pre_gCell){
  fCellScoreComposer.EstimatorCalculation(aStep);
  ScorePopulation(pre_gCell, aStep);
}

void G4CellScorer::ScoreAnEnteringStep(const G4Step &aStep,
				       const G4GeometryCell &post_gCell){
  // population counting
  ScorePopulation(post_gCell, aStep);
  // entering tracks
  fCellScoreComposer.TrackEnters();
}

void G4CellScorer::ScoreAnInVolumeStep(const G4Step &aStep,
				       const G4GeometryCell &post_gCell){
  // population counting
  ScorePopulation(post_gCell, aStep);
  // collisions counting
  if (aStep.GetPostStepPoint()->GetStepStatus() != fGeomBoundary) {
    fCellScoreComposer.SetCollisionWeight(aStep.GetTrack()->GetWeight());
  }

  // counting for step length estimators
  fCellScoreComposer.EstimatorCalculation(aStep);

}

void G4CellScorer::ScorePopulation(const G4GeometryCell &post_gCell, 
				   const G4Step &aStep) {
  //    check for new event
  fTrackLogger.SetEventID(G4EventManager::GetEventManager()->
			  GetConstCurrentEvent()->
			  GetEventID());
  //    increase population
  if (fTrackLogger.FirstEnterance(aStep.GetTrack()->GetTrackID())) {
    fCellScoreComposer.NewTrackPopedUp();
  }
}

const G4CellScoreComposer &G4CellScorer::GetCellScoreComposer() const {
  return fCellScoreComposer;
}

const G4CellScoreValues &G4CellScorer::GetCellScoreValues() const {
  return fCellScoreComposer.GetStandardCellScoreValues();
}

