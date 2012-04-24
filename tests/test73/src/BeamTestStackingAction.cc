//
//  BeamTestStackingAction.cc
//  MSCTest
//
//  Created by Andrea Dotti on 3/9/12.
//  Copyright (c) 2012 CERN. All rights reserved.
//


#include <iostream>
#include "BeamTestStackingAction.hh"
#include "G4Track.hh"
#include "G4EventManager.hh"
#include "BeamTestEventAction.hh"

BeamTestStackingAction::BeamTestStackingAction() : 
    evtmgr(G4EventManager::GetEventManager())
{
    
}

G4ClassificationOfNewTrack BeamTestStackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
    G4ClassificationOfNewTrack result( fUrgent );
    if ( aTrack->GetParentID() != 0 ) {
        //This is a secondary... Kill it and kill events,
        //a secondary has been produced...
        //G4cout<<"Secondary ("<<aTrack->GetParentID()<<") Creator";
        //const G4VProcess* proc = aTrack->GetCreatorProcess();
        //if ( proc ) G4cout<<proc->GetProcessName()<<G4endl;
        //else G4cout<<"Unknown process"<<G4endl;
        result = fKill;
        evtmgr->AbortCurrentEvent();
        //Signal to EventAction that event is aborted, hits collected
        //so far will not be procssed...
        G4UserEventAction* theEventAction=const_cast<G4UserEventAction*>(evtmgr->GetUserEventAction());
        static_cast<BeamTestEventAction*>(theEventAction)->AbortEvent();
        //G4cout<<"AborEventRequested"<<G4endl;
    }
    return result;
}