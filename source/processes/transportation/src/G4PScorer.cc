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
// $Id: G4PScorer.cc,v 1.7 2002-04-09 17:40:15 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4PScorer.cc
//
// ----------------------------------------------------------------------

#include "G4PScorer.hh"
#include "G4Step.hh"
#include "G4PStep.hh"
#include "G4Sigma.hh"
#include "G4StateManager.hh"

G4PScorer::G4PScorer(){}
G4PScorer::~G4PScorer(){}

void G4PScorer::Score(const G4Step &aStep, const G4PStep &aPstep)
{
  G4Track *track = aStep.GetTrack();

  if (track->GetTrackStatus()==fStopAndKill) {
    G4cout << "      track status is StopAndKill -> do nothing" << G4endl;
  }
  // do the scoring
  else {    
    G4double weight = track->GetWeight();
    // the map fPtkTallys, the map nametallys and the "G4Sigma" will
    // be setup the first time they are accesed. 
    G4PTouchableKey post_ptk(aPstep.fPostTouchableKey); 
    if (aPstep.fCrossBoundary) { 
      // Pstep crosses boundary
      fPtkTallys[post_ptk]["WeighteOfHistorysEntering"].Xin(weight);
      fPtkTallys[post_ptk]["HistorysEnteringWeighted"].Xin(1,weight);
      fPtkTallys[post_ptk]["EnergyEnteringHistoryWeighted"].
	Xin(track->GetKineticEnergy(), weight);
    } 
    else { 
      // Pstep with both points in the same I volume
      if (aStep.GetPostStepPoint()->GetStepStatus() != fGeomBoundary) {
	fPtkTallys[post_ptk]["WeighteOfCollisions"].Xin(weight);
	fPtkTallys[post_ptk]["CollisionsWeighted"].Xin(1,weight);
      }
    }
  }
}

G4std::ostream& operator<<(G4std::ostream &out, const G4PScorer &ps)
{
  out << ps.GetMapPtkTallys();
  return out;
}
