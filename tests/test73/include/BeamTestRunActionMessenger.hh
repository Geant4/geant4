//
//  BeamTestRunActionMessenger.hh
//  MSCTest
//
//  Created by Andrea Dotti on 3/13/12.
//  Copyright (c) 2012 CERN. All rights reserved.
//

#ifndef MSCTest_BeamTestRunActionMessenger_hh
#define MSCTest_BeamTestRunActionMessenger_hh

#include "G4UImessenger.hh"

class BeamTestRunAction;
class G4UIdirectory;
class G4UIcmdWithAString;

class BeamTestRunActionMessenger : public G4UImessenger {
public:
    BeamTestRunActionMessenger( BeamTestRunAction* );
    virtual ~BeamTestRunActionMessenger();
    void SetNewValue( G4UIcommand*, G4String );
private:
    BeamTestRunAction* theRunAction;
    G4UIdirectory* theDir;
    G4UIcmdWithAString* theSuffixCmd;
    G4UIcmdWithAString* theMergedFileCmd;
    G4UIcmdWithAString* theMergeCmd;
    G4UIcmdWithAString* theFinalizeCmd;
};

#endif
