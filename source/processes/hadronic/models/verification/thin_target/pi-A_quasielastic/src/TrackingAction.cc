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
// $Id: TrackingAction.cc,v 1.1 2003-07-31 01:22:14 dwright Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "TrackingAction.hh"
#include "EventAction.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"


TrackingAction::TrackingAction(EventAction* theAction)
  : evtAction(theAction)
{ }


TrackingAction::~TrackingAction()
{ }


void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Store each secondary produced by primary track, then kill secondary track

  if (aTrack->GetParentID() == 1) {
    evtAction->StoreDynamicParticle(aTrack->GetDynamicParticle() );
    fpTrackingManager->GetTrack()->SetTrackStatus(fStopAndKill);
  }
}
