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
// $Id: G4CellScorer.cc,v 1.4 2006-06-29 18:09:31 gunter Exp $
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

void G4CellScorer::ScorePopulation(const G4GeometryCell &, 
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

