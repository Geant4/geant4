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
// $Id: G4PIScorer.cc,v 1.4 2002-04-09 17:40:15 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4PIScorer.cc
//
// ----------------------------------------------------------------------

#include "G4PIScorer.hh"
#include "G4Step.hh"
#include "G4PStep.hh"
#include "G4Sigma.hh"
#include "G4StateManager.hh"
#include "G4VIStore.hh"
#include "G4PStepStream.hh"
#include "G4StepPoint.hh"

G4PIScorer::G4PIScorer(const G4VIStore &IStore)
  : fIStore(IStore), fNerr(0),fCorrectWeight(true)
{}

G4PIScorer::~G4PIScorer()
{}

void G4PIScorer::Score(const G4Step &aStep, const G4PStep &aPstep)
{
  G4Track *track = aStep.GetTrack();

  if (track->GetTrackStatus()==fStopAndKill) {
    G4cout << "      track status is StopAndKill -> do nothing" << G4endl;
  }
  // do the scoring
  else {    
    G4double weight = track->GetWeight();
    G4double import = fIStore.GetImportance(aPstep.fPostTouchableKey);
    G4double i_w =import  * weight;
    if (i_w != 1.) {
      fNerr++;
      G4StepPoint* prepoint = aStep.GetPreStepPoint();
      G4StepPoint* postpoint = aStep.GetPostStepPoint();
      
      G4cout << "G4PIScorer::Score:" << G4endl;
      G4cout << "  Stepnumber = " << aStep.GetTrack()->GetCurrentStepNumber() 
	     << G4endl;
      G4cout << "  TrackID = " << aStep.GetTrack()->GetTrackID() << G4endl;
      G4cout << "  ParentID = " << aStep.GetTrack()->GetParentID() << G4endl;
      G4cout << "  Track wieght = " << weight << G4endl;
      G4cout << "  Importance = " << import << G4endl;
      G4cout << "  i_w = " << i_w << G4endl;
      G4cout << "  aPstep: " << aPstep << G4endl;
      G4cout << "  steplength = " << aStep.GetStepLength() << G4endl;
      G4cout << "  pre_pos = " << prepoint->GetPosition()
	     << ", pre_dir = " << prepoint->GetMomentumDirection() << G4endl;
      G4cout << "  post_pos = " << postpoint->GetPosition()
	     << ", post_dir = " << postpoint->GetMomentumDirection() << G4endl;
      G4cout << "  vxt mom: " << track->GetVertexMomentumDirection() << G4endl;
      fCorrectWeight = false;
    }
    fPScorer.Score(aStep, aPstep);
  }
}

G4std::ostream& operator<<(G4std::ostream &out, const G4PIScorer &ps)
{
  out << ps.GetMapPtkTallys();
  return out;
}
