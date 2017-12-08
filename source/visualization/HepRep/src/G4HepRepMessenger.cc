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
// $Id: G4HepRepMessenger.cc 106384 2017-10-09 09:34:51Z gcosmo $
//
#include "G4HepRepMessenger.hh"

G4HepRepMessenger*
G4HepRepMessenger::fpInstance = 0;

G4HepRepMessenger*
G4HepRepMessenger::GetInstance()
{
	if (!fpInstance) fpInstance = new G4HepRepMessenger;
	return fpInstance;
}

G4HepRepMessenger::G4HepRepMessenger() :
    fileDir(""),
	fileName("G4Data"),
	overwrite(false),
	cullInvisibles(false),
	cylAsPolygons(false),
	scale(1.),
    suffix (""),
    geometry(true),
    pointAttributes(false),
    solids(true),
    invisibles(true) {

    heprepDirectory = new G4UIdirectory("/vis/heprep/");
    heprepDirectory->SetGuidance("HepRep commands.");
		
	setFileDirCommand = new G4UIcmdWithAString("/vis/heprep/setFileDir", this);
	setFileDirCommand->SetGuidance("Set directory for output.");
	setFileDirCommand->SetGuidance("This command is used by HepRepFile, not by HepRepXML.");
	setFileDirCommand->SetParameterName("directory",false);
	if ( getenv( "G4HEPREPFILE_DIR" ) == NULL ) {
		setFileDirCommand->SetDefaultValue("");
	} else {
		setFileDirCommand->SetDefaultValue(getenv("G4HEPREPFILE_DIR"));
		fileDir = getenv("G4HEPREPFILE_DIR");
	}
	setFileDirCommand->AvailableForStates(G4State_Idle);
		
	setFileNameCommand = new G4UIcmdWithAString("/vis/heprep/setFileName", this);
	setFileNameCommand->SetGuidance("Set file name for output.");
	setFileNameCommand->SetGuidance("This command is used by HepRepFile, not by HepRepXML.");
	setFileNameCommand->SetParameterName("directory",false);
	if ( getenv( "G4HEPREPFILE_NAME" ) == NULL ) {
		setFileNameCommand->SetDefaultValue("G4Data");
	} else {
		setFileNameCommand->SetDefaultValue(getenv("G4HEPREPFILE_NAME"));
		fileName = getenv("G4HEPREPFILE_NAME");
	}
	setFileNameCommand->AvailableForStates(G4State_Idle);

	setOverwriteCommand = new G4UIcmdWithABool("/vis/heprep/setOverwrite", this);
	setOverwriteCommand->SetGuidance("Set true to write all output to exact same file name.");
	setOverwriteCommand->SetGuidance("Set false to increment the file name for each new output.");
	setOverwriteCommand->SetGuidance("This command is used by HepRepFile, not by HepRepXML.");
	setOverwriteCommand->SetParameterName("flag",false);
	if ( getenv( "G4HEPREPFILE_OVERWRITE" ) == NULL ) {
		setOverwriteCommand->SetDefaultValue(false);
	} else {
		setOverwriteCommand->SetDefaultValue(getenv("G4HEPREPFILE_OVERWRITE"));
		overwrite = setOverwriteCommand->ConvertToBool(getenv("G4HEPREPFILE_OVERWRITE"));
	}
	setOverwriteCommand->AvailableForStates(G4State_Idle);
		
	setCullInvisiblesCommand = new G4UIcmdWithABool("/vis/heprep/setCullInvisibles", this);
	setCullInvisiblesCommand->SetGuidance("Remove invisible objects from output file.");
	setCullInvisiblesCommand->SetGuidance("This command is used by HepRepFile, not by HepRepXML.");
	setCullInvisiblesCommand->SetParameterName("flag",false);
	if ( getenv( "G4HEPREPFILE_CULL" ) == NULL ) {
		setCullInvisiblesCommand->SetDefaultValue(false);
	} else {
		setCullInvisiblesCommand->SetDefaultValue(getenv("G4HEPREPFILE_CULL"));
		cullInvisibles = setCullInvisiblesCommand->ConvertToBool(getenv("G4HEPREPFILE_CULL"));
	}
	setCullInvisiblesCommand->AvailableForStates(G4State_Idle);
		
	renderCylAsPolygonsCommand = new G4UIcmdWithABool("/vis/heprep/renderCylAsPolygons", this);
	renderCylAsPolygonsCommand->SetGuidance("Render cylinders and cones as polygons.");
	renderCylAsPolygonsCommand->SetGuidance("This command is used by HepRepFile, not by HepRepXML.");
	renderCylAsPolygonsCommand->SetParameterName("flag",false);
	renderCylAsPolygonsCommand->SetDefaultValue(false);
	renderCylAsPolygonsCommand->AvailableForStates(G4State_Idle);
		
	setScaleCommand = new G4UIcmdWithADouble("/vis/heprep/scale",this);
	setScaleCommand->SetGuidance("Re-Scale coordinates.");
	setScaleCommand->SetParameterName("Scale",true);
	setScaleCommand->SetDefaultValue(1.);
	setScaleCommand->SetRange("Scale > 0");
	
	setCenterCommand = new G4UIcmdWith3VectorAndUnit("/vis/heprep/center",this);
	setCenterCommand->SetGuidance("Re-Center coordinates.");
	setCenterCommand->SetParameterName("CenterX","CenterY","CenterZ",true);
	setCenterCommand->SetDefaultValue(G4ThreeVector(0.,0.,0.));
	setCenterCommand->SetDefaultUnit("m");
		
    setEventNumberSuffixCommand = new G4UIcmdWithAString("/vis/heprep/setEventNumberSuffix", this);
    setEventNumberSuffixCommand->SetGuidance("Write separate event files, appended with given suffix.");
    setEventNumberSuffixCommand->SetGuidance("Define the suffix with a pattern such as '-0000'.");
	setEventNumberSuffixCommand->SetGuidance("This command is used by HepRepXML, not by HepRepFile.");
    setEventNumberSuffixCommand->SetParameterName("suffix",false);
    setEventNumberSuffixCommand->SetDefaultValue("");
    setEventNumberSuffixCommand->AvailableForStates(G4State_Idle);
    
    appendGeometryCommand = new G4UIcmdWithABool("/vis/heprep/appendGeometry", this);
    appendGeometryCommand->SetGuidance("Appends copy of geometry to every event.");
	appendGeometryCommand->SetGuidance("This command is used by HepRepXML, not by HepRepFile.");
    appendGeometryCommand->SetParameterName("flag",false);
    appendGeometryCommand->SetDefaultValue(true);
    appendGeometryCommand->AvailableForStates(G4State_Idle);

    addPointAttributesCommand = new G4UIcmdWithABool("/vis/heprep/addPointAttributes", this);
    addPointAttributesCommand->SetGuidance("Adds point attributes to the points of trajectories.");
	addPointAttributesCommand->SetGuidance("This command is used by HepRepXML, not by HepRepFile.");
    addPointAttributesCommand->SetParameterName("flag",false);
    addPointAttributesCommand->SetDefaultValue(false);
    addPointAttributesCommand->AvailableForStates(G4State_Idle);
		
	useSolidsCommand = new G4UIcmdWithABool("/vis/heprep/useSolids", this);
	useSolidsCommand->SetGuidance("Use HepRep Solids, rather than Geant4 Primitives.");
	useSolidsCommand->SetGuidance("This command is used by HepRepXML, not by HepRepFile..");
	useSolidsCommand->SetParameterName("flag",false);
	useSolidsCommand->SetDefaultValue(true);
	useSolidsCommand->AvailableForStates(G4State_Idle);
}

G4HepRepMessenger::~G4HepRepMessenger() {
	delete setFileDirCommand;
	delete setFileNameCommand;
	delete setOverwriteCommand;
	delete setCullInvisiblesCommand;
    delete renderCylAsPolygonsCommand;
	delete setScaleCommand;
	delete setCenterCommand;
    delete setEventNumberSuffixCommand;
    delete appendGeometryCommand;
    delete addPointAttributesCommand;
    delete useSolidsCommand;
    delete heprepDirectory;
}

G4String G4HepRepMessenger::GetCurrentValue(G4UIcommand * command) {
    if (command==setFileDirCommand) {
        return fileDir;
    } else if (command==setFileNameCommand) {
        return fileName; 
    } else if (command==setOverwriteCommand) {
        return overwrite; 
    } else if (command==setCullInvisiblesCommand) {
        return cullInvisibles; 
    } else if (command==renderCylAsPolygonsCommand) {
        return renderCylAsPolygonsCommand->ConvertToString(cylAsPolygons);
    } else if (command==setScaleCommand) {
        return setScaleCommand->ConvertToString(scale);
    } else if (command==setCenterCommand) {
        return setCenterCommand->ConvertToString(center,"m");
    } else if (command==setEventNumberSuffixCommand) {
        return suffix; 
    } else if (command==appendGeometryCommand) {
        return appendGeometryCommand->ConvertToString(geometry); 
    } else if (command==addPointAttributesCommand) {
        return addPointAttributesCommand->ConvertToString(pointAttributes); 
    } else if (command==useSolidsCommand) {
        return useSolidsCommand->ConvertToString(solids);
    } else {
        return "";
    }
}

void G4HepRepMessenger::SetNewValue(G4UIcommand * command, G4String newValue) {
    if (command==setFileDirCommand) {
        fileDir = newValue;
    } else if (command==setFileNameCommand) {
        fileName = newValue;
    } else if (command==setOverwriteCommand) {
        overwrite = setOverwriteCommand->GetNewBoolValue(newValue);
    } else if (command==setCullInvisiblesCommand) {
		cullInvisibles = setCullInvisiblesCommand->GetNewBoolValue(newValue);
    } else if (command==renderCylAsPolygonsCommand) {
        cylAsPolygons = renderCylAsPolygonsCommand->GetNewBoolValue(newValue);
    } else if (command==setScaleCommand) {
        scale = setScaleCommand->GetNewDoubleValue(newValue);
    } else if (command==setCenterCommand) {
        center = setCenterCommand->GetNew3VectorValue(newValue);
    } else if (command==setEventNumberSuffixCommand) {
        suffix = newValue;
    } else if (command==appendGeometryCommand) {
        geometry = appendGeometryCommand->GetNewBoolValue(newValue);
    } else if (command==addPointAttributesCommand) {
        pointAttributes = addPointAttributesCommand->GetNewBoolValue(newValue);
    } else if (command==useSolidsCommand) {
        solids = useSolidsCommand->GetNewBoolValue(newValue);
    } 
}

G4String G4HepRepMessenger::getFileDir() {
    return fileDir;
}

G4String G4HepRepMessenger::getFileName() {
    return fileName;
}

G4bool G4HepRepMessenger::getOverwrite() {
    return overwrite;
}

G4bool G4HepRepMessenger::getCullInvisibles() {
    return cullInvisibles;
}

G4bool G4HepRepMessenger::renderCylAsPolygons() {
    return cylAsPolygons;
}

G4double G4HepRepMessenger::getScale() {
    return scale;
}

G4ThreeVector G4HepRepMessenger::getCenter() {
    return center;
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

G4bool G4HepRepMessenger::useSolids() {
    return solids;
}

G4bool G4HepRepMessenger::writeInvisibles() {
    return invisibles;
}

