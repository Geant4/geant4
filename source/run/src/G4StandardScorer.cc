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
// $Id: G4StandardScorer.cc,v 1.1 2002-07-11 16:19:44 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4StandardScorer.hh"

#include "G4Step.hh"
#include "G4PStep.hh"
#include "G4StateManager.hh"
#include "G4PStepStream.hh"
#include "G4VIStore.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

G4StandardScorer::G4StandardScorer()
{;}

G4StandardScorer::~G4StandardScorer()
{;}

void G4StandardScorer::Score(const G4Step &aStep, const G4PStep &aPstep)
{
  G4Track *track = aStep.GetTrack();
  
  if (track->GetTrackStatus()==fStopAndKill) {
    G4cout << "      track status is StopAndKill -> do nothing" << G4endl;
  }
  else {    // do the scoring
    G4PTouchableKey pre_ptk(aPstep.fPreTouchableKey); 
    G4PTouchableKey post_ptk(aPstep.fPostTouchableKey); 
    G4PTouchableKey *sl_ptk = 0;

    fMapPtkTrackLogger[post_ptk].
      SetEventID(G4EventManager::GetEventManager()->
		 GetConstCurrentEvent()->
		 GetEventID());
    
    // create the raw values for the step length estimators
    // for a cell
    G4SLRawValues slr(aStep);
    
    // find which cell to score 
    // and score non sl scores
    if (aPstep.fCrossBoundary) { 
      // the sl scores belong to the pre_ptk cell
      sl_ptk = &pre_ptk;
      // track enters the post_ptk cell
      fPtkScores[post_ptk].TrackEnters();
      // check if track entered the post_ptk cell
      // and increment population if track is new
      if (fMapPtkTrackLogger[post_ptk].
	  FirstEnterance(track->GetTrackID())) {
	fPtkScores[post_ptk].NewTrackPopedUp();
      }
    } 
    else { 
      // the sl scores belong to the post_cell
      sl_ptk = &post_ptk;
      // add up tracks created inside post ptk cell to population
      if  (track->GetCurrentStepNumber() == 1) {
	fPtkScores[post_ptk].NewTrackPopedUp();
      }
      
      // score collisions 
      // don't count collisions with mass geometry bounderies
      if (aStep.GetPostStepPoint()->GetStepStatus() != fGeomBoundary) {
	fPtkScores[post_ptk].SetCollisionWeight(track->GetWeight());
      }
      
    }
    // score the sl estimators
    fPtkScores[*sl_ptk].SetSLRawValues(slr);
  }
  

  
}

G4std::ostream& operator<<(G4std::ostream &out, 
			   const G4StandardScorer &ps)
{
  out << ps.GetMapPtkStandardCellScorer();
  return out;
}


