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
// $Id: SteppingAction.cc,v 1.1 2003-07-31 01:21:51 dwright Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4Track.hh"


SteppingAction::SteppingAction(EventAction* theAction)
  : evtAction(theAction)
{ }


SteppingAction::~SteppingAction()
{ }


void SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  G4Track* aTrack = aStep->GetTrack();

  // For primary track only
  if (aTrack->GetTrackID() == 1) {

    // If momentum changes in step, store particle, kill track
 
    G4double zdiff = abs(aStep->GetPreStepPoint()->GetMomentum().z() - 
                         aStep->GetPostStepPoint()->GetMomentum().z());
 
    if (zdiff > 1.e-10) {
      evtAction->StoreDynamicParticle(aTrack->GetDynamicParticle());
      aTrack->SetTrackStatus(fStopAndKill);
      G4cout << " Primary track, ID = "
             << aTrack->GetTrackID() << " stored " << G4endl;
    }
  }
}


