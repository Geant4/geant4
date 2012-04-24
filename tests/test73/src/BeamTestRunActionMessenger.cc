//
//  BeamTestRunActionMessenger.cc
//  MSCTest
//
//  Created by Andrea Dotti on 3/13/12.
//  Copyright (c) 2012 CERN. All rights reserved.
//

#include <iostream>
#include "BeamTestRunAction.hh"
#include "BeamTestRunActionMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

BeamTestRunActionMessenger::BeamTestRunActionMessenger( BeamTestRunAction* runAction ) : theRunAction(runAction)
{
    theDir = new G4UIdirectory( "/analysis/" );
    theDir->SetGuidance("Commands controlling the analysis");
    
    theSuffixCmd = new G4UIcmdWithAString("/analysis/filename/addSuffix",this);
    theSuffixCmd->SetGuidance("Add a suffix to the output filename.");
    theSuffixCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    
    theMergedFileCmd = new G4UIcmdWithAString("/analysis/filename/mergedFileName",this);
    theMergedFileCmd->SetGuidance("File name for the merged file name (see /analysis/filename/merge command).");
    theMergedFileCmd->SetDefaultValue("merged.root");
    
    theMergeCmd = new G4UIcmdWithAString("/analysis/mergeFiles",this);
    theMergeCmd->SetGuidance("Merge multiple file names matching the pattern in a single file");
    theMergeCmd->SetDefaultValue("*.root");
    
    theFinalizeCmd = new G4UIcmdWithAString("/analysis/rootmacro",this);
    theFinalizeCmd->SetGuidance("Execute ROOT macro.");
    theFinalizeCmd->SetDefaultValue("macro.C");
    
}

BeamTestRunActionMessenger::~BeamTestRunActionMessenger() 
{
    delete theDir;
    delete theSuffixCmd;
    delete theMergedFileCmd;
    delete theMergeCmd;
    delete theFinalizeCmd;
}

void BeamTestRunActionMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
    if ( cmd == theSuffixCmd ) {
        theRunAction->AddSuffix( newValue );
    }
    if ( cmd == theMergedFileCmd ) {
        theRunAction->SetMergedFilename(newValue);
    }
    if ( cmd == theMergeCmd ) {
        theRunAction->MergeFiles(newValue);
    }
    if ( cmd == theFinalizeCmd ) {
        theRunAction->Finalize(newValue);
    }
}