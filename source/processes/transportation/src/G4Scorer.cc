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
// $Id: G4Scorer.cc,v 1.1 2002-06-13 07:12:33 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4Scorer.hh"

#include "G4Step.hh"
#include "G4PStep.hh"
#include "G4Sigma.hh"
#include "G4StateManager.hh"
#include "G4PStepStream.hh"
#include "G4VIStore.hh"

G4Scorer::G4Scorer()
{;}

G4Scorer::~G4Scorer()
{;}

void G4Scorer::Score(const G4Step &aStep, const G4PStep &aPstep)
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
    // the map fPtkTallys, the map nametallys and the "G4Sigma" will
    // be setup the first time they are accesed. 
    G4PTouchableKey pre_ptk(aPstep.fPreTouchableKey); 
    G4PTouchableKey post_ptk(aPstep.fPostTouchableKey); 
    
    if (aPstep.fCrossBoundary)
    { 
      // Pstep crosses boundary
      fPtkTallys[post_ptk]["WeighteOfHistorysEntering"].Xin(weight);
      fPtkTallys[post_ptk]["HistorysEnteringWeighted"].Xin(1,weight);
      fPtkTallys[post_ptk]["EnergyEnteringHistoryWeighted"].
	Xin(track->GetKineticEnergy(), weight);
      // add step length in of previous cell 
      fPtkTallys[pre_ptk]["StepLentgh"].Xin(aStep.GetStepLength(),weight);
      // add energy times step length. 
      fPtkTallys[pre_ptk]["SteplengthTimesEnergy"]
	.Xin(track->GetKineticEnergy()*aStep.GetStepLength(),
	     weight);
    } 
    else
      { 
	// add step length to the current volume
	fPtkTallys[post_ptk]["StepLentgh"].
	  Xin(aStep.GetStepLength(),weight);
	// add energy times weight. 
	fPtkTallys[post_ptk]["SteplengthTimesEnergy"]
	  .Xin(track->GetKineticEnergy()*aStep.GetStepLength(),
	       weight);
	// Pstep with both points in the same I volume
	if (aStep.GetPostStepPoint()->GetStepStatus() != fGeomBoundary)
	  {
	    fPtkTallys[post_ptk]["WeighteOfCollisions"].Xin(weight);
	    fPtkTallys[post_ptk]["CollisionsWeighted"].Xin(1,weight);
	    fPtkTallys[post_ptk]["CollisionEnergyWeighted"].
	      Xin(track->GetKineticEnergy(),weight);
	  }
      }
  }
}

G4std::ostream& operator<<(G4std::ostream &out, const G4Scorer &ps)
{
  out << ps.GetMapPtkTallys();
  return out;
}
