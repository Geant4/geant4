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
//  BeamTestRunActionMessenger.cc
//  MSCTest
//
//  Created by Andrea Dotti on 3/13/12.
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
