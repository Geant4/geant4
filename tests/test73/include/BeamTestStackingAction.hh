//
//  BeamTestStackingAction.hh
//  MSCTest
//
//  Created by Andrea Dotti on 3/9/12.
//  Copyright (c) 2012 CERN. All rights reserved.
//

#ifndef MSCTest_BeamTestStackingAction_hh
#define MSCTest_BeamTestStackingAction_hh

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ClassificationOfNewTrack.hh"

class G4Track;
class G4EventManager;

class BeamTestStackingAction : public G4UserStackingAction {
public:
    BeamTestStackingAction();
    virtual ~BeamTestStackingAction() {};
    virtual G4ClassificationOfNewTrack ClassifyNewTrack( const G4Track* aTrack);
private:
    G4EventManager* evtmgr;
};

#endif
