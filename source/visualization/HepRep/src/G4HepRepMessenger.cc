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
#include "G4HepRepMessenger.hh"

G4HepRepMessenger::G4HepRepMessenger() :
    suffix (""),
    geometry(true) {

    heprepDirectory = new G4UIdirectory("/vis/heprep/");
    heprepDirectory->SetGuidance("HepRep commands.");

    setEventNumberSuffixCommand = new G4UIcmdWithAString("/vis/heprep/setEventNumberSuffix", this);
    setEventNumberSuffixCommand->SetGuidance("Write separate event files, appended with given suffix.");
    setEventNumberSuffixCommand->SetGuidance("Define the suffix with a pattern such as '-0000'.");
    setEventNumberSuffixCommand->SetParameterName("suffix",false);
    setEventNumberSuffixCommand->SetDefaultValue("");
    setEventNumberSuffixCommand->AvailableForStates(G4State_Idle);
    
/*
    appendGeometryCommand = new G4UIcmdWithABool("/vis/heprep/appendGeometry", this);
    appendGeometryCommand->SetGuidance("Appends copy of geometry to every event.");
    appendGeometryCommand->SetParameterName("flag",false);
    appendGeometryCommand->SetDefaultValue(true);
    appendGeometryCommand->AvailableForStates(G4State_Idle);
*/

    addPointAttributesCommand = new G4UIcmdWithABool("/vis/heprep/addPointAttributes", this);
    addPointAttributesCommand->SetGuidance("Adds point attributes to the points of trajectories.");
    addPointAttributesCommand->SetParameterName("flag",false);
    addPointAttributesCommand->SetDefaultValue(false);
    addPointAttributesCommand->AvailableForStates(G4State_Idle);
}

G4HepRepMessenger::~G4HepRepMessenger() {
    delete setEventNumberSuffixCommand;
//    delete appendGeometryCommand;
    delete addPointAttributesCommand;
    delete heprepDirectory;
}

G4String G4HepRepMessenger::GetCurrentValue(G4UIcommand * command) {
    if (command==setEventNumberSuffixCommand) {
        return suffix;
//    } else if (command==appendGeometryCommand) {
//        return appendGeometryCommand->ConvertToString(geometry); 
    } else if (command==addPointAttributesCommand) {
        return addPointAttributesCommand->ConvertToString(pointAttributes); 
    } else {
        return "";
    }
}

void G4HepRepMessenger::SetNewValue(G4UIcommand * command, G4String newValue) {
    if (command==setEventNumberSuffixCommand) {
        suffix = newValue;
//    } else if (command==appendGeometryCommand) {
//        geometry = appendGeometryCommand->GetNewBoolValue(newValue);
    } else if (command==addPointAttributesCommand) {
        pointAttributes = addPointAttributesCommand->GetNewBoolValue(newValue);
    }
}

G4String G4HepRepMessenger::getEventNumberSuffix() {
    return suffix;
}

G4bool G4HepRepMessenger::appendGeometry() {
    return geometry;
}

G4bool G4HepRepMessenger::addPointAttributes() {
    return pointAttributes;
}
