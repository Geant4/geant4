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
// $Id: B02Scorer.cc,v 1.3 2002-04-19 10:54:27 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "B02Scorer.hh"
#include "G4Step.hh"
#include "G4PStep.hh"
#include "G4Sigma.hh"
#include "G4StateManager.hh"

B02Scorer::B02Scorer(){}
B02Scorer::~B02Scorer(){}

void B02Scorer::Score(const G4Step &aStep, const G4PStep &aPstep)
{
  G4Track *track = aStep.GetTrack();

  if (track->GetTrackStatus()==fStopAndKill)
  {
    G4cout << "      track status is StopAndKill -> do nothing" << G4endl;
  }
  // do the scoring
  else
  {    
    // the map fPtkTallys, the map nametallys and the "G4Sigma" will
    // be setup the first time they are accesed. 
    // the user may choos any other place to store the results e.g. 
    // histogramms
    G4PTouchableKey post_ptk(aPstep.fPostTouchableKey); 
    if (aPstep.fCrossBoundary)
    { 
      // Pstep crosses boundary
      fPtkTallys[post_ptk]["HistorysEntering"].Xin(1);
      fPtkTallys[post_ptk]["EnergyEnteringHistory"].
	Xin(track->GetKineticEnergy());
    } 
    else
    { 
      // this check is necessary in case of scoring in a parallel world
      if (aStep.GetPostStepPoint()->GetStepStatus() != fGeomBoundary)
      {
	fPtkTallys[post_ptk]["Collisions"].Xin(1);
      }
    }
  }
}

G4std::ostream& operator<<(G4std::ostream &out, const B02Scorer &ps)
{
  out << ps.GetMapPtkTallys();
  return out;
}
