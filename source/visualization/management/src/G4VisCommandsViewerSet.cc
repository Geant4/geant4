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
// $Id: G4VisCommandsViewerSet.cc,v 1.17 2001-11-12 18:22:12 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/viewer/set commands - John Allison  16th May 2000

#include "G4VisCommandsViewerSet.hh"

#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UnitsTable.hh"
#include "G4VisManager.hh"

G4VisCommandsViewerSet::G4VisCommandsViewerSet ():
  fLightsVector    (G4ThreeVector(1.,1.,1.)),
  fUpVector        (G4ThreeVector(0.,1.,0.)),
  fViewpointVector (G4ThreeVector(0.,0.,1.))
{
  G4bool omitable;
  G4UIparameter* parameter;

  fpCommandAll = new G4UIcmdWithAString ("/vis/viewer/set/all",this);
  fpCommandAll->SetGuidance
    ("/vis/viewer/set/all <from-viewer-name>"
     "\nCopies view parameters from from-viewer to current viewer.");
  fpCommandAll->SetParameterName ("from-viewer-name",omitable = false);
  viewerNameCommands.push_back (fpCommandAll);

  fpCommandAutoRefresh = new G4UIcmdWithABool
    ("/vis/viewer/set/autoRefresh",this);
  fpCommandAutoRefresh->SetGuidance
    ("/vis/viewer/set/autoRefresh [true|false]");
  fpCommandAutoRefresh->SetGuidance
    ("View is automatically refreshed after a change of view parameters.");
  fpCommandAutoRefresh->SetParameterName("auto-refresh",omitable = true);
  fpCommandAutoRefresh->SetDefaultValue(false);

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
  fpCommandLightsMove->SetParameterName("lightsMove",omitable = false);
  // fpCommandLightsMove->SetCandidates("move-with-camera move-with-object");
  // Own parsing.

  fpCommandLightsThetaPhi = new G4UIcommand
    ("/vis/viewer/set/lightsThetaPhi", this);
  fpCommandLightsThetaPhi -> SetGuidance
    ("/vis/viewer/set/lightsThetaPhi  [<theta>] [<phi>] [deg|rad]");
  fpCommandLightsThetaPhi -> SetGuidance
    ("Set direction from target to lights.");
  parameter = new G4UIparameter("theta", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandLightsThetaPhi -> SetParameter (parameter);
  parameter = new G4UIparameter("phi", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandLightsThetaPhi -> SetParameter (parameter);
  parameter = new G4UIparameter ("unit", 's', omitable = true);
  parameter -> SetDefaultValue ("deg");
  fpCommandLightsThetaPhi -> SetParameter (parameter);

  fpCommandLightsVector = new G4UIcommand
    ("/vis/viewer/set/lightsVector", this);
  fpCommandLightsVector -> SetGuidance
    ("/vis/viewer/set/lightsVector  [<x>] [<y>] [<z>]");
  fpCommandLightsVector -> SetGuidance
    ("Set direction from target to lights.");
  parameter = new G4UIparameter("x", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandLightsVector -> SetParameter (parameter);
  parameter = new G4UIparameter("y", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandLightsVector -> SetParameter (parameter);
  parameter = new G4UIparameter ("z", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandLightsVector -> SetParameter (parameter);

  fpCommandLineSegments = new G4UIcmdWithAnInteger
    ("/vis/viewer/set/lineSegmentsPerCircle",this);
  fpCommandLineSegments->SetGuidance
    ("/vis/viewer/set/lineSegmentsPerCircle  [<number-of-sides-per-circle>]");
  fpCommandLineSegments->SetGuidance
    ("  Number of sides per circle in polygon/polyhedron graphical"
     "\nrepresentation of objects with curved lines/surfaces.");
  fpCommandLineSegments->SetParameterName("line-segments",omitable = true);
  fpCommandLineSegments->SetDefaultValue(24);

  fpCommandProjection = new G4UIcommand("/vis/viewer/set/projection",this);
  fpCommandProjection->SetGuidance
    ("/vis/viewer/set/projection"
     " o[rthogonal]|p[erspective] [<field-half-angle>] [deg|rad]");
  fpCommandProjection->SetGuidance
    ("Note: 1st parameter will be parsed for first character \"o\" or \"p\".");
  fpCommandProjection->SetGuidance("Default: orthogonal 30 deg");
  parameter = new G4UIparameter("projection",'s',omitable = true);
  parameter->SetDefaultValue("orthogonal");
  fpCommandProjection->SetParameter(parameter);
  parameter = new G4UIparameter("field-half-angle",'d',omitable = true);
  parameter->SetDefaultValue(30.);
  fpCommandProjection->SetParameter(parameter);
  parameter = new G4UIparameter("unit",'s',omitable = true);
  parameter->SetDefaultValue("deg");
  fpCommandProjection->SetParameter(parameter);

  fpCommandSectionPlane = new G4UIcommand 
    ("/vis/viewer/set/sectionPlane",this);
  fpCommandSectionPlane -> SetGuidance
    (
     "Set plane for drawing section (DCUT).  Specify plane by"
     "\nx y z units nx ny nz, e.g., for a y-z plane at x = 1 cm:"
     "\n/vis/viewer/set/sectionPlane on 1 0 0 cm 1 0 0"
     );
  parameter  =  new G4UIparameter("Selector",'c',true);
  parameter  -> SetDefaultValue  ("?");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("x",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Coordinate of point on the plane.");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("y",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Coordinate of point on the plane.");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("z",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Coordinate of point on the plane.");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("unit",'s',omitable = true);
  parameter  -> SetDefaultValue  ("cm");
  parameter  -> SetGuidance      ("Unit of point on the plane.");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("nx",'d',omitable = true);
  parameter  -> SetDefaultValue  (1);
  parameter  -> SetGuidance      ("Component of plane normal.");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("ny",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Component of plane normal.");
  fpCommandSectionPlane->SetParameter(parameter);
  parameter  =  new G4UIparameter("nz",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Component of plane normal.");
  fpCommandSectionPlane->SetParameter(parameter);

  fpCommandStyle = new G4UIcmdWithAString ("/vis/viewer/set/style",this);
  fpCommandStyle->SetGuidance ("/vis/viewer/set/style wireframe|surface");
  fpCommandStyle->SetGuidance
    ("Note: parameter will be parsed for first character \"w\" or \"s\".");
  fpCommandStyle->SetParameterName ("style",omitable = false);
  // fpCommandStyle->SetCandidates("wireframe surface");  // Own parsing.

  fpCommandUpThetaPhi = new G4UIcommand
    ("/vis/viewer/set/upThetaPhi", this);
  fpCommandUpThetaPhi -> SetGuidance
    ("/vis/viewer/set/upThetaPhi  [<theta>] [<phi>] [deg|rad]");
  fpCommandUpThetaPhi -> SetGuidance
    ("Set up vector.  Viewer will attempt always to show"
     " this direction upwards.");
  parameter = new G4UIparameter("theta", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandUpThetaPhi -> SetParameter (parameter);
  parameter = new G4UIparameter("phi", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandUpThetaPhi -> SetParameter (parameter);
  parameter = new G4UIparameter ("unit", 's', omitable = true);
  parameter -> SetDefaultValue ("deg");
  fpCommandUpThetaPhi -> SetParameter (parameter);

  fpCommandUpVector = new G4UIcommand
    ("/vis/viewer/set/upVector", this);
  fpCommandUpVector -> SetGuidance
    ("/vis/viewer/set/upVector  [<x>] [<y>] [<z>]");
  fpCommandUpVector -> SetGuidance
    ("Set up vector.  Viewer will attempt always to show"
     " this direction upwards.");
  parameter = new G4UIparameter("x", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandUpVector -> SetParameter (parameter);
  parameter = new G4UIparameter("y", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandUpVector -> SetParameter (parameter);
  parameter = new G4UIparameter ("z", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandUpVector -> SetParameter (parameter);

  fpCommandViewpointThetaPhi = new G4UIcommand
    ("/vis/viewer/set/viewpointThetaPhi", this);
  fpCommandViewpointThetaPhi -> SetGuidance
    ("/vis/viewer/set/viewpointThetaPhi  [<theta>] [<phi>] [deg|rad]");
  fpCommandViewpointThetaPhi -> SetGuidance
    ("Set direction from target to camera.  Also changes lightpoint direction"
     "\nif lights are set to move with camera.");
  parameter = new G4UIparameter("theta", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandViewpointThetaPhi -> SetParameter (parameter);
  parameter = new G4UIparameter("phi", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandViewpointThetaPhi -> SetParameter (parameter);
  parameter = new G4UIparameter ("unit", 's', omitable = true);
  parameter -> SetDefaultValue ("deg");
  fpCommandViewpointThetaPhi -> SetParameter (parameter);

  fpCommandViewpointVector = new G4UIcommand
    ("/vis/viewer/set/viewpointVector", this);
  fpCommandViewpointVector -> SetGuidance
    ("/vis/viewer/set/viewpointVector  [<x>] [<y>] [<z>]");
  fpCommandViewpointVector -> SetGuidance
    ("Set direction from target to camera.  Also changes lightpoint direction"
     "\nif lights are set to move with camera.");
  parameter = new G4UIparameter("x", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandViewpointVector -> SetParameter (parameter);
  parameter = new G4UIparameter("y", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandViewpointVector -> SetParameter (parameter);
  parameter = new G4UIparameter ("z", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandViewpointVector -> SetParameter (parameter);
}

G4VisCommandsViewerSet::~G4VisCommandsViewerSet() {
  delete fpCommandAll;
  delete fpCommandAutoRefresh;
  delete fpCommandCulling;
  delete fpCommandEdge;
  delete fpCommandHiddenEdge;
  delete fpCommandHiddenMarker;
  delete fpCommandLineSegments;
  delete fpCommandLightsMove;
  delete fpCommandLightsThetaPhi;
  delete fpCommandLightsVector;
  delete fpCommandProjection;
  delete fpCommandSectionPlane;
  delete fpCommandStyle;
  delete fpCommandUpThetaPhi;
  delete fpCommandUpVector;
  delete fpCommandViewpointThetaPhi;
  delete fpCommandViewpointVector;
}

G4String G4VisCommandsViewerSet::GetCurrentValue(G4UIcommand* command) {
  G4String currentValue;
  if (command == fpCommandLightsThetaPhi) {
    currentValue = ConvertToString(fLightsVector.theta(),
			   fLightsVector.phi(), "deg");
  }
  else if (command == fpCommandLightsVector) {
    currentValue = ConvertToString(fLightsVector);
  }
  else if (command == fpCommandViewpointThetaPhi) {
    currentValue = ConvertToString(fViewpointVector.theta(),
			   fViewpointVector.phi(), "deg");
  }
  else if (command == fpCommandViewpointVector) {
    currentValue = ConvertToString(fViewpointVector);
  }
  else {
    currentValue = "invalid";
  }
  return currentValue;
}

void G4VisCommandsViewerSet::SetNewValue
(G4UIcommand* command,G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* currentViewer = fpVisManager->GetCurrentViewer();
  if (!currentViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cout << 
	"ERROR: G4VisCommandsViewerSet::SetNewValue: no current viewer."
	     << G4endl;
    }
    return;
  }

  G4ViewParameters vp = currentViewer->GetViewParameters();

  if (command == fpCommandAll) {
    G4VViewer* fromViewer = fpVisManager->GetViewer(newValue);
    if (!fromViewer) {
      if (verbosity >= G4VisManager::errors) {
	G4cout <<
	  "ERROR: G4VisCommandsViewerSet::SetNewValue: all:"
	  "\n  unrecognised from-viewer."
	       << G4endl;
      }
      return;
    }
    if (fromViewer == currentViewer) {
      if (verbosity >= G4VisManager::warnings) {
      G4cout <<
	"WARNING: G4VisCommandsViewerSet::SetNewValue: all:"
	"\n  from-viewer and current viewer are identical."
	     << G4endl;
      }
      return;
    }
    vp = fromViewer->GetViewParameters();
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "View parameters of viewer \"" << currentViewer->GetName()
	     << "\"\n  set to those of viewer \"" << fromViewer->GetName()
	     << "\"."
	     << G4endl;
    }
  }

  else if (command == fpCommandAutoRefresh) {
    G4bool autoRefresh = GetNewBoolValue(newValue);
    vp.SetAutoRefresh(autoRefresh);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Views will ";
      if (!vp.IsAutoRefresh()) G4cout << "not ";
      G4cout << "be automatically refreshed after a change of view parameters."
	     << G4endl;
    }
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
      if (verbosity >= G4VisManager::confirmations) {
	G4cout <<
	  "G4VisCommandsViewerSet::SetNewValue: culling: global culling flag"
	  " set to " << ConvertToString(boolFlag) <<
	  ".\n  Does not change specific culling flags."
	       << G4endl;
      }
    }
    else if (cullingOption == "coveredDaughters") {
      vp.SetCullingCovered(boolFlag);
      if (verbosity >= G4VisManager::confirmations) {
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
    }
    else if (cullingOption == "invisible") {
      vp.SetCullingInvisible(boolFlag);
      if (verbosity >= G4VisManager::confirmations) {
	G4cout <<
	  "G4VisCommandsViewerSet::SetNewValue: culling: culling invisible"
	  "\n  flag set to "
	       << boolFlag << ConvertToString(boolFlag) <<
	  ".  Volumes marked invisible will be culled,"
	  "\n  i.e., not drawn, if this flag is true."
	       << G4endl;
      }
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
      if (verbosity >= G4VisManager::confirmations) {
	G4cout <<
	  "G4VisCommandsViewerSet::SetNewValue: culling: culling by density"
	  "\n  flag set to " << ConvertToString(boolFlag) <<
	  ".  Volumes with density less than " <<
	  G4BestUnit(density,"Volumic Mass") <<
	  "\n  will be culled, i.e., not drawn, if this flag is true."
	       << G4endl;
      }
    }
    else {
      if (verbosity >= G4VisManager::errors) {
	G4cout <<
	  "ERROR: G4VisCommandsViewerSet::SetNewValue: culling:"
	  "\n  option not recognised."
	       << G4endl;
      }
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
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Drawing style of viewer \"" << currentViewer->GetName()
	     << "\" set to " << vp.GetDrawingStyle()
	     << G4endl;
    }
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
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Drawing style of viewer \"" << currentViewer->GetName()
	     << "\" set to " << vp.GetDrawingStyle()
	     << G4endl;
    }
  }

  else if (command == fpCommandHiddenMarker) {
    G4bool hidden = GetNewBoolValue(newValue);
    if (hidden) vp.SetMarkerHidden();
    else vp.SetMarkerNotHidden();
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Markers will ";
      if (vp.IsMarkerNotHidden()) G4cout << "not ";
      G4cout << "be hidden under solid objects." << G4endl;
    }
  }

  else if (command == fpCommandLightsMove) {
    G4String s (newValue);
    if (s.find("cam") != G4String::npos) vp.SetLightsMoveWithCamera(true);
    else if(s.find("obj") != G4String::npos) vp.SetLightsMoveWithCamera(false);
    else {
      if (verbosity >= G4VisManager::errors) {
	G4cout << "ERROR: \"" << newValue << "\" not recognised."
	"  Looking for \"cam\" or \"obj\" in string." << G4endl;
      }
    }
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Lights move with ";
      if (vp.GetLightsMoveWithCamera())
	G4cout << "camera (object appears to rotate).";
      else G4cout << "object (the viewer appears to be moving).";
      G4cout << G4endl;
    }
  }

  else if (command == fpCommandLightsThetaPhi) {
    G4double theta, phi;
    GetNewDoublePairValue(newValue, theta, phi);
    G4double x = sin (theta) * cos (phi);
    G4double y = sin (theta) * sin (phi);
    G4double z = cos (theta);
    fLightsVector = G4ThreeVector (x, y, z);
    vp.SetLightpointDirection(fLightsVector);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Lights direction set to "
	     << vp.GetLightpointDirection() << G4endl;
    }
  }

  else if (command == fpCommandLightsVector) {
    fLightsVector = GetNew3VectorValue(newValue);
    vp.SetLightpointDirection(fLightsVector);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Lights direction set to "
	     << vp.GetLightpointDirection() << G4endl;
    }
  }

  else if (command == fpCommandLineSegments) {
    G4int nSides = GetNewIntValue(newValue);
    vp.SetNoOfSides(nSides);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout <<
	"Number of line segements per circle in polygon approximation is "
	     << nSides << G4endl;
    }
  }

  else if (command == fpCommandProjection) {
    G4double fieldHalfAngle;
    size_t iPos = 0;
    if (newValue[iPos] == 'o') {  // "orthogonal"
      fieldHalfAngle = 0.;
    }
    else if (newValue[iPos] == 'p') {  // "perspective"
      G4String dummy;
      G4String unit;
      G4std::istrstream is ((char*)newValue.data());
      is >> dummy >> fieldHalfAngle >> unit;
      fieldHalfAngle *= ValueOf(unit);
      if (fieldHalfAngle > 89.5 * deg || fieldHalfAngle <= 0.0) {
	if (verbosity >= G4VisManager::errors) {
	  G4cout <<
	    "ERROR: Field half angle should be 0 < angle <= 89.5 degrees.";
	  G4cout << G4endl;
	}
	return;
      }
    }
    else {
      if (verbosity >= G4VisManager::errors) {
	G4cout << "ERROR: \"" << newValue << "\" not recognised."
	  "  Looking for 'o' or 'p' first character." << G4endl;
      }
      return;
    }
    vp.SetFieldHalfAngle(fieldHalfAngle);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Projection style of viewer \"" << currentViewer->GetName()
	     << "\" set to ";
      if (fieldHalfAngle == 0.) {
	G4cout << "orthogonal.";
      }
      else {
	G4cout << "perspective\n  with half angle " << fieldHalfAngle / deg
	       << " degrees.";
      }
      G4cout << G4endl;
    }
  }

  else if (command == fpCommandSectionPlane) {
    G4String choice, unit;
    G4double x, y, z, nx, ny, nz;
    const char* t = newValue;
    G4std::istrstream is ((char*)t);
    is >> choice >> x >> y >> z >> unit >> nx >> ny >> nz;

    G4int iSelector = -1;
    if (choice.compareTo("off",G4String::ignoreCase) == 0) iSelector = 0;
    if (choice.compareTo("on",G4String::ignoreCase) == 0) iSelector = 1;
    if (iSelector < 0) {
      if (verbosity >= G4VisManager::errors) {
	G4cout << "Choice not recognised (on/off)." << G4endl;
	G4cout << "Section drawing is currently: ";
	if (vp.IsSection ()) G4cout << "on";
	else                    G4cout << "off";
	G4cout << ".\nSection plane is currently: "
	       << vp.GetSectionPlane ();
	G4cout << G4endl;
      }
      return;
    }

    G4double F;
    switch (iSelector) {
    case 0:
      vp.UnsetSectionPlane();
      break;
    case 1:
      F = ValueOf(unit);
      x *= F; y *= F; z *= F;
      vp.SetSectionPlane(G4Plane3D(G4Normal3D(nx,ny,nz),
				   G4Point3D(x,y,z)));
      vp.SetViewpointDirection(G4Normal3D(nx,ny,nz));
      break;
    default:
      if (verbosity >= G4VisManager::errors) {
	G4cout << "ERROR: Choice not recognised (on/off)."
	       << G4endl;
      }
      break;
    }

    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Section drawing is now: ";
      if (vp.IsSection ()) G4cout << "on";
      else                    G4cout << "off";
      G4cout << ".\nSection plane is now: "
	     << vp.GetSectionPlane ();
      G4cout << G4endl;
    }
  }

  else if (command == fpCommandStyle) {
    G4ViewParameters::DrawingStyle existingStyle = vp.GetDrawingStyle();
    size_t iPos = 0;
    if (newValue[iPos] == 'w') {  // "wireframe"
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
    else if (newValue[iPos] == 's') {  // "surface"
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
      if (verbosity >= G4VisManager::errors) {
	G4cout << "ERROR: \"" << newValue << "\" not recognised."
	  "  Looking for 'w' or 's' first character." << G4endl;
      }
      return;
    }
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Drawing style of viewer \"" << currentViewer->GetName()
	     << "\" set to " << vp.GetDrawingStyle()
	     << G4endl;
    }
  }

  else if (command == fpCommandUpThetaPhi) {
    G4double theta, phi;
    GetNewDoublePairValue(newValue, theta, phi);
    G4double x = sin (theta) * cos (phi);
    G4double y = sin (theta) * sin (phi);
    G4double z = cos (theta);
    fUpVector = G4ThreeVector (x, y, z);
    vp.SetUpVector(fUpVector);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Up direction set to " << vp.GetUpVector() << G4endl;
    }
  }

  else if (command == fpCommandUpVector) {
    fUpVector = GetNew3VectorValue(newValue);
    vp.SetUpVector(fUpVector);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Up direction set to " << vp.GetUpVector() << G4endl;
    }
  }

  else if (command == fpCommandViewpointThetaPhi) {
    G4double theta, phi;
    GetNewDoublePairValue(newValue, theta, phi);
    G4double x = sin (theta) * cos (phi);
    G4double y = sin (theta) * sin (phi);
    G4double z = cos (theta);
    fViewpointVector = G4ThreeVector (x, y, z);
    vp.SetViewAndLights(fViewpointVector);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Viewpoint direction set to "
	     << vp.GetViewpointDirection() << G4endl;
      if (vp.GetLightsMoveWithCamera ()) {
	G4cout << "Lightpoint direction set to "
	       << vp.GetActualLightpointDirection () << G4endl;
      }
    }
  }

  else if (command == fpCommandViewpointVector) {
    fViewpointVector = GetNew3VectorValue(newValue);
    vp.SetViewAndLights(fViewpointVector);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Viewpoint direction set to "
	     << vp.GetViewpointDirection() << G4endl;
      if (vp.GetLightsMoveWithCamera ()) {
	G4cout << "Lightpoint direction set to "
	       << vp.GetActualLightpointDirection () << G4endl;
      }
    }
  }

  else {
    if (verbosity >= G4VisManager::errors) {
      G4cout <<
	"ERROR: G4VisCommandsViewerSet::SetNewValue: unrecognised command."
	     << G4endl;
    }
    return;
  }

  SetViewParameters(currentViewer,vp);
}
