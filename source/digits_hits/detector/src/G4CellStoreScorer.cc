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
// $Id: G4CellStoreScorer.cc,v 1.2 2006/06/29 18:05:43 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4CellStoreScorer
//
// ----------------------------------------------------------------------

#include "G4CellStoreScorer.hh"

#include "G4VCellScorer.hh"

#include "G4Track.hh"
#include "G4GeometryCell.hh"
#include "G4Step.hh"
#include "G4GeometryCellStep.hh"
#include "G4VCellScorerStore.hh"

G4CellStoreScorer::G4CellStoreScorer(G4VCellScorerStore &csc) :
  fCellScorerStore(csc)
{}

G4CellStoreScorer::~G4CellStoreScorer()
{}


void G4CellStoreScorer::
Score(const G4Step &aStep, const G4GeometryCellStep &aPstep){
  G4Track *track = 0;
  track = aStep.GetTrack();
  
  if (track->GetTrackStatus()==fStopAndKill) {
    G4cout << "G4CellStoreScorer::Score:   track status is StopAndKill -> do nothing" << G4endl;
  }
  else { 
    // chose the cells to be scored
    G4GeometryCell pre_gCell(aPstep.GetPreGeometryCell()); 
    G4GeometryCell post_gCell(aPstep.GetPostGeometryCell()); 
    G4VCellScorer *post_cs = 0;
    post_cs = fCellScorerStore.GetCellScore(post_gCell);
    if (aPstep.GetCrossBoundary()) { 
      // entering post_gCell
      if (post_cs) {
	post_cs->ScoreAnEnteringStep(aStep, post_gCell);
      }
      // exiting pre_gCell
      G4VCellScorer *pre_cs = 0;
      pre_cs = fCellScorerStore.GetCellScore(pre_gCell);
      if (pre_cs) {
	pre_cs->ScoreAnExitingStep(aStep, pre_gCell);
      }
    } 
    else {
      // step in post_gCell
      if (post_cs) {
	post_cs->ScoreAnInVolumeStep(aStep, post_gCell);
      }
    }
  }
  
}

