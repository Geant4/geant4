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
// $Id: B06Scorer.cc,v 1.3 2002/04/19 10:54:33 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

#include "B06Scorer.hh"

#include "G4Step.hh"
#include "G4PStep.hh"
#include "G4Sigma.hh"
#include "G4StateManager.hh"
#include "G4PStepStream.hh"
#include "G4VIStore.hh"

B06Scorer::B06Scorer(const G4VIStore &IStore)
  : fIStore(IStore), fCorrectWeight(true)
{;}

B06Scorer::~B06Scorer()
{;}

void B06Scorer::Score(const G4Step &aStep, const G4PStep &aPstep)
{
  G4Track *track = aStep.GetTrack();

  if (track->GetTrackStatus()==fStopAndKill)
  {
    G4cout << "      track status is StopAndKill -> do nothing" << G4endl;
  }
  // do the scoring
  else
  {    
    G4double weight = track->GetWeight();
    // check if importances are "correct" -------------------------------
    G4double import = fIStore.GetImportance(aPstep.fPostTouchableKey);
    G4double i_w =import  * weight;
    if (i_w != 1.) {
      G4StepPoint* prepoint = aStep.GetPreStepPoint();
      G4StepPoint* postpoint = aStep.GetPostStepPoint();
      
      G4cout << "B06Scorer::Score:" << G4endl;
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
    // --------------------------------------------------------------------

    
    // the map fPtkTallys, the map nametallys and the "G4Sigma" will
    // be setup the first time they are accesed. 
    G4PTouchableKey post_ptk(aPstep.fPostTouchableKey); 
    if (aPstep.fCrossBoundary)
    { 
      // Pstep crosses boundary
      fPtkTallys[post_ptk]["WeighteOfHistorysEntering"].Xin(weight);
      fPtkTallys[post_ptk]["HistorysEnteringWeighted"].Xin(1,weight);
      fPtkTallys[post_ptk]["EnergyEnteringHistoryWeighted"].
	Xin(track->GetKineticEnergy(), weight);
    } 
    else
    { 
      // Pstep with both points in the same I volume
      if (aStep.GetPostStepPoint()->GetStepStatus() != fGeomBoundary)
      {
	fPtkTallys[post_ptk]["WeighteOfCollisions"].Xin(weight);
	fPtkTallys[post_ptk]["CollisionsWeighted"].Xin(1,weight);
      }
    }
  }
}

G4std::ostream& operator<<(G4std::ostream &out, const B06Scorer &ps)
{
  out << ps.GetMapPtkTallys();
  return out;
}
