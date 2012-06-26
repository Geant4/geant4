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
//  BeamTestStackingAction.cc
//  MSCTest
//
//  Created by Andrea Dotti on 3/9/12.
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
