// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsViewerSet.cc,v 1.4 2001-02-04 01:37:38 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/viewer/set commands - John Allison  16th May 2000

#include "G4VisCommandsViewerSet.hh"

#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UnitsTable.hh"
#include "G4VisManager.hh"

G4VisCommandsViewerSet::G4VisCommandsViewerSet () {
  G4bool omitable;
  G4UIparameter* parameter;

  fpCommandAll = new G4UIcmdWithAString ("/vis/viewer/set/all",this);
  fpCommandAll->SetGuidance
    ("/vis/viewer/set/all <from-viewer-name>"
     "\nCopies view parameters from from-viewer to current viewer.");
  fpCommandAll->SetParameterName ("from-viewer-name",omitable = false);
  viewerNameCommands.push_back (fpCommandAll);

  fpCommandCulling = new G4UIcommand("/vis/viewer/set/culling",this);
  fpCommandCulling->SetGuidance
    ("/vis/viewer/set/culling global|coveredDaughters|invisible|density"
     " [true|false] [density] [unit]");
  fpCommandCulling->SetGuidance
    ("If density option is true, provide threshold density and unit");
  fpCommandCulling->SetGuidance
    ("Unit is g/cm3, mg/cm3 or kg/m3");
  parameter = new G4UIparameter("culling-option",'s',omitable = false);
  parameter->SetParameterCandidates
    ("global coveredDaughters invisible density");
  fpCommandCulling->SetParameter(parameter);
  parameter = new G4UIparameter("flag",'b',omitable = true);
  parameter->SetDefaultValue("true");
  fpCommandCulling->SetParameter(parameter);
  parameter = new G4UIparameter("density",'d',omitable = true);
  parameter->SetDefaultValue("0.01");
  fpCommandCulling->SetParameter(parameter);
  parameter = new G4UIparameter("unit",'s',omitable = true);
  parameter->SetDefaultValue("g/cm3");
  fpCommandCulling->SetParameter(parameter);

  fpCommandEdge = new G4UIcmdWithABool("/vis/viewer/set/edge",this);
  fpCommandEdge->SetGuidance("/vis/viewer/set/edge [true|false]");
  fpCommandEdge->SetParameterName("edge",omitable = true);
  fpCommandEdge->SetDefaultValue(true);

  fpCommandHiddenEdge =
    new G4UIcmdWithABool("/vis/viewer/set/hiddenEdge",this);
  fpCommandHiddenEdge->SetGuidance("/vis/viewer/set/edge [true|false]");
  fpCommandHiddenEdge->SetParameterName("hidden-edge",omitable = true);
  fpCommandHiddenEdge->SetDefaultValue(true);

  fpCommandHiddenMarker =
    new G4UIcmdWithABool("/vis/viewer/set/hiddenMarker",this);
  fpCommandHiddenMarker->SetGuidance
    ("/vis/viewer/set/hiddenMarker [true|false]");
  fpCommandHiddenMarker->SetParameterName("hidden-marker",omitable = true);
  fpCommandHiddenMarker->SetDefaultValue(true);

  fpCommandLightsMove = new G4UIcmdWithAString
    ("/vis/viewer/set/lightsMove",this);
  fpCommandLightsMove->SetGuidance
    ("/vis/viewer/set/lightsMove with-camera|with-object");
  fpCommandLightsMove->SetGuidance
    ("Note: parameter will be parsed for \"cam\" or \"obj\".");
  fpCommandLightsMove->SetParameterName ("lightsMove",omitable = false);
  // fpCommandLightsMove->SetCandidates("move-with-camera move-with-object");
  // Own parsing.

  fpCommandStyle = new G4UIcmdWithAString ("/vis/viewer/set/style",this);
  fpCommandStyle->SetGuidance ("/vis/viewer/set/style wireframe|surface");
  fpCommandStyle->SetGuidance
    ("Note: parameter will be parsed for first character \"w\" or \"s\".");
  fpCommandStyle->SetParameterName ("style",omitable = false);
  // fpCommandStyle->SetCandidates("wireframe surface");  // Own parsing.
}

G4VisCommandsViewerSet::~G4VisCommandsViewerSet() {
  delete fpCommandAll;
  delete fpCommandCulling;
  delete fpCommandEdge;
  delete fpCommandHiddenEdge;
  delete fpCommandHiddenMarker;
  delete fpCommandStyle;
}

G4String G4VisCommandsViewerSet::GetCurrentValue(G4UIcommand* command) {
  return "invalid";
}

void G4VisCommandsViewerSet::SetNewValue
(G4UIcommand* command,G4String newValue) {

  G4VViewer* currentViewer = fpVisManager->GetCurrentViewer();
  if (!currentViewer) {
    G4cout << "G4VisCommandsViewerSet::SetNewValue: no current viewer."
	   << G4endl;
    return;
  }

  G4ViewParameters vp = currentViewer->GetViewParameters();

  if (command == fpCommandAll) {
    G4VViewer* fromViewer = fpVisManager->GetViewer(newValue);
    if (!fromViewer) {
      G4cout <<
	"G4VisCommandsViewerSet::SetNewValue: all: unrecognised from-viewer."
	     << G4endl;
      return;
    }
    if (fromViewer == currentViewer) {
      G4cout <<
	"G4VisCommandsViewerSet::SetNewValue: all:"
	"\n  from-viewer and current viewer are identical."
	     << G4endl;
      return;
    }
    vp = fromViewer->GetViewParameters();
    G4cout << "View parameters of viewer \"" << currentViewer->GetName()
	   << "\"\n  set to those of viewer \"" << fromViewer->GetName()
	   << "\"."
	   << G4endl;
  }

  else if (command == fpCommandCulling) {
    G4String cullingOption, stringFlag, unit;
    G4double density;
    G4bool boolFlag;
    G4std::istrstream is ((char*)newValue.data());
    is >> cullingOption >> stringFlag >> density >> unit;
    boolFlag = GetNewBoolValue(stringFlag);
    if (cullingOption == "global") {
      vp.SetCulling(boolFlag);
      G4cout <<
	"G4VisCommandsViewerSet::SetNewValue: culling: global culling flag"
	" set to " << ConvertToString(boolFlag) <<
	".\n  Does not change specific culling flags."
	     << G4endl;
    }
    else if (cullingOption == "coveredDaughters") {
      vp.SetCullingCovered(boolFlag);
      G4cout <<
	"G4VisCommandsViewerSet::SetNewValue: culling: culling covered"
	"\n  daughters flag set to "
	     << ConvertToString(boolFlag) <<
	".  Daughters covered by opaque mothers"
	"\n  will be culled, i.e., not drawn, if this flag is true."
	"\n  Note: this is only effective in surface drawing style,"
	"\n  and then only if the volumes are visible and opaque, and then"
	"\n  only if no sections or cutways are in operation."
	     << G4endl;
    }
    else if (cullingOption == "invisible") {
      vp.SetCullingInvisible(boolFlag);
      G4cout <<
	"G4VisCommandsViewerSet::SetNewValue: culling: culling invisible"
	"\n  flag set to "
	     << boolFlag << ConvertToString(boolFlag) <<
	".  Volumes marked invisible will be culled,"
	"\n  i.e., not drawn, if this flag is true."
	     << G4endl;
    }
    else if (cullingOption == "density") {
      vp.SetDensityCulling(boolFlag);
      if (boolFlag) {
	density *= G4UnitDefinition::GetValueOf(unit);
	vp.SetVisibleDensity(density);
      }
      else {
	density = vp.GetVisibleDensity();
      }
      G4cout <<
	"G4VisCommandsViewerSet::SetNewValue: culling: culling by density"
	"\n  flag set to " << ConvertToString(boolFlag) <<
	".  Volumes with density less than " <<
	G4BestUnit(density,"Volumic Mass") <<
	"\n  will be culled, i.e., not drawn, if this flag is true."
	     << G4endl;
    }
    else {
      G4cout <<
	"G4VisCommandsViewerSet::SetNewValue: culling: option not recognised."
	     << G4endl;
    }
  }

  else if (command == fpCommandEdge) {
    G4ViewParameters::DrawingStyle existingStyle = vp.GetDrawingStyle();
    if (GetNewBoolValue(newValue)) {
      switch (existingStyle) {
      case G4ViewParameters::wireframe:
	break;
      case G4ViewParameters::hlr:
	break;
      case G4ViewParameters::hsr:
	vp.SetDrawingStyle(G4ViewParameters::hlhsr);
        break;
      case G4ViewParameters::hlhsr:
        break;
      }
    }
    else {
      switch (existingStyle) {
      case G4ViewParameters::wireframe:
	break;
      case G4ViewParameters::hlr:
	break;
      case G4ViewParameters::hsr:
	break;
      case G4ViewParameters::hlhsr:
	vp.SetDrawingStyle(G4ViewParameters::hsr);
	break;
      }
    }
    G4cout << "Drawing style of viewer \"" << currentViewer->GetName()
	   << "\" set to " << vp.GetDrawingStyle()
	   << G4endl;
  }

  else if (command == fpCommandHiddenEdge) {
    G4ViewParameters::DrawingStyle existingStyle = vp.GetDrawingStyle();
    if (GetNewBoolValue(newValue)) {
      switch (existingStyle) {
      case G4ViewParameters::wireframe:
	vp.SetDrawingStyle(G4ViewParameters::hlr);
	break;
      case G4ViewParameters::hlr:
	break;
      case G4ViewParameters::hsr:
	vp.SetDrawingStyle(G4ViewParameters::hlhsr);
        break;
      case G4ViewParameters::hlhsr:
        break;
      }
    }
    else {
      switch (existingStyle) {
      case G4ViewParameters::wireframe:
	break;
      case G4ViewParameters::hlr:
	vp.SetDrawingStyle(G4ViewParameters::wireframe);
	break;
      case G4ViewParameters::hsr:
	break;
      case G4ViewParameters::hlhsr:
	vp.SetDrawingStyle(G4ViewParameters::hsr);
	break;
      }
    }
    G4cout << "Drawing style of viewer \"" << currentViewer->GetName()
	   << "\" set to " << vp.GetDrawingStyle()
	   << G4endl;
  }

  else if (command == fpCommandHiddenMarker) {
    G4bool hidden = GetNewBoolValue(newValue);
    if (hidden) vp.SetMarkerHidden();
    else vp.SetMarkerNotHidden();
    G4cout << "Markers will ";
    if (vp.IsMarkerNotHidden()) G4cout << "not ";
    G4cout << "be hidden under solid objects." << G4endl;
  }

  else if (command == fpCommandStyle) {
    G4ViewParameters::DrawingStyle existingStyle = vp.GetDrawingStyle();
    if (newValue[0] == 'w') {  // "wireframe"
      switch (existingStyle) {
      case G4ViewParameters::wireframe:
	break;
      case G4ViewParameters::hlr:
	break;
      case G4ViewParameters::hsr:
	vp.SetDrawingStyle(G4ViewParameters::wireframe);
        break;
      case G4ViewParameters::hlhsr:
	vp.SetDrawingStyle(G4ViewParameters::hlr);
        break;
      }
    }
    else if (newValue[0] == 's') {  // "surface"
      switch (existingStyle) {
      case G4ViewParameters::wireframe:
	vp.SetDrawingStyle(G4ViewParameters::hsr);
	break;
      case G4ViewParameters::hlr:
	vp.SetDrawingStyle(G4ViewParameters::hlhsr);
	break;
      case G4ViewParameters::hsr:
	break;
      case G4ViewParameters::hlhsr:
	break;
      }
    }
    else {
      G4cout <<
	"G4VisCommandsViewerSet::SetNewValue: style: unrecognised style."
	     << G4endl;
      return;
    }
    G4cout << "Drawing style of viewer \"" << currentViewer->GetName()
	   << "\" set to " << vp.GetDrawingStyle()
	   << G4endl;
  }

  else if (command == fpCommandLightsMove) {
    G4String s (newValue);
    if (s.find("cam") != G4String::npos) vp.SetLightsMoveWithCamera(true);
    else if(s.find("obj") != G4String::npos) vp.SetLightsMoveWithCamera(false);
    else {
      G4cout << "\"" << newValue << "\" not recognised."
	"  Looking for \"cam\" or \"obj\" in string." << G4endl;
    }
    G4cout << "Lights move with ";
    if (vp.GetLightsMoveWithCamera())
      G4cout << "camera (object appears to rotate).";
    else G4cout << "object (the viewer appears to be moving).";
    G4cout << G4endl;
  }

  else {
    G4cout << "G4VisCommandsViewerSet::SetNewValue: unrecognised command."
	   << G4endl;
    return;
  }

  currentViewer->SetViewParameters(vp);
  G4cout << "Issue /vis/viewer/refresh to see effect." << G4endl;
  // For now...
  fpVisManager->SetCurrentViewParameters() = vp;
}
