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
// $Id: G4VisCommandsViewer.cc 104163 2017-05-15 06:52:42Z gcosmo $

// /vis/viewer commands - John Allison  25th October 1998

#include "G4VisCommandsViewer.hh"

#include "G4VisManager.hh"
#include "G4GraphicsSystemList.hh"
#include "G4VisCommandsScene.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4Point3D.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#ifdef G4VIS_USE_STD11
#include <chrono>
#include <thread>
#endif
#include <sstream>
#include <fstream>
#include <sstream>
#include <iomanip>
#ifdef WIN32
#include <regex>
#include <filesystem>
#endif //WIN32

G4VVisCommandViewer::G4VVisCommandViewer () {}

G4VVisCommandViewer::~G4VVisCommandViewer () {}

void G4VVisCommandViewer::SetViewParameters
(G4VViewer* viewer, const G4ViewParameters& viewParams) {
  viewer->SetViewParameters(viewParams);
  RefreshIfRequired(viewer);
}

void G4VVisCommandViewer::RefreshIfRequired(G4VViewer* viewer) {
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4VSceneHandler* sceneHandler = viewer->GetSceneHandler();
  const G4ViewParameters& viewParams = viewer->GetViewParameters();
  if (sceneHandler && sceneHandler->GetScene()) {
    if (viewParams.IsAutoRefresh()) {
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
    }
    else {
      if (verbosity >= G4VisManager::warnings) {
	G4cout << "Issue /vis/viewer/refresh or flush to see effect."
	       << G4endl;
      }
    }
  }
}

////////////// /vis/viewer/addCutawayPlane ///////////////////////////////////////

G4VisCommandViewerAddCutawayPlane::G4VisCommandViewerAddCutawayPlane () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/viewer/addCutawayPlane", this);
  fpCommand -> SetGuidance
    ("Add cutaway plane to current viewer.");
  G4UIparameter* parameter;
  parameter  =  new G4UIparameter("x",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Coordinate of point on the plane.");
  fpCommand->SetParameter(parameter);
  parameter  =  new G4UIparameter("y",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Coordinate of point on the plane.");
  fpCommand->SetParameter(parameter);
  parameter  =  new G4UIparameter("z",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Coordinate of point on the plane.");
  fpCommand->SetParameter(parameter);
  parameter  =  new G4UIparameter("unit",'s',omitable = true);
  parameter  -> SetDefaultValue  ("m");
  parameter  -> SetGuidance      ("Unit of point on the plane.");
  fpCommand->SetParameter(parameter);
  parameter  =  new G4UIparameter("nx",'d',omitable = true);
  parameter  -> SetDefaultValue  (1);
  parameter  -> SetGuidance      ("Component of plane normal.");
  fpCommand->SetParameter(parameter);
  parameter  =  new G4UIparameter("ny",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Component of plane normal.");
  fpCommand->SetParameter(parameter);
  parameter  =  new G4UIparameter("nz",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Component of plane normal.");
  fpCommand->SetParameter(parameter);
}

G4VisCommandViewerAddCutawayPlane::~G4VisCommandViewerAddCutawayPlane () {
  delete fpCommand;
}

G4String G4VisCommandViewerAddCutawayPlane::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandViewerAddCutawayPlane::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
  "ERROR: No current viewer - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }

  G4double x, y, z, nx, ny, nz;
  G4String unit;
  std::istringstream is (newValue);
  is >> x >> y >> z >> unit >> nx >> ny >> nz;
  G4double F = G4UIcommand::ValueOf(unit);
  x *= F; y *= F; z *= F;

  G4ViewParameters vp = viewer->GetViewParameters();
  vp.AddCutawayPlane(G4Plane3D(G4Normal3D(nx,ny,nz), G4Point3D(x,y,z)));
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Cutaway planes for viewer \"" << viewer->GetName() << "\" now:";
    const G4Planes& cutaways = vp.GetCutawayPlanes();
    for (size_t i = 0; i < cutaways.size(); ++i)
      G4cout << "\n  " << i << ": " << cutaways[i];
    G4cout << G4endl;
  }

  SetViewParameters(viewer, vp);
}

////////////// /vis/viewer/changeCutawayPlane ///////////////////////////////////////

G4VisCommandViewerChangeCutawayPlane::G4VisCommandViewerChangeCutawayPlane () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/viewer/changeCutawayPlane", this);
  fpCommand -> SetGuidance("Change cutaway plane.");
  G4UIparameter* parameter;
  parameter  =  new G4UIparameter("index",'i',omitable = false);
  parameter  -> SetGuidance      ("Index of plane: 0, 1, 2.");
  fpCommand->SetParameter(parameter);
  parameter  =  new G4UIparameter("x",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Coordinate of point on the plane.");
  fpCommand->SetParameter(parameter);
  parameter  =  new G4UIparameter("y",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Coordinate of point on the plane.");
  fpCommand->SetParameter(parameter);
  parameter  =  new G4UIparameter("z",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Coordinate of point on the plane.");
  fpCommand->SetParameter(parameter);
  parameter  =  new G4UIparameter("unit",'s',omitable = true);
  parameter  -> SetDefaultValue  ("m");
  parameter  -> SetGuidance      ("Unit of point on the plane.");
  fpCommand->SetParameter(parameter);
  parameter  =  new G4UIparameter("nx",'d',omitable = true);
  parameter  -> SetDefaultValue  (1);
  parameter  -> SetGuidance      ("Component of plane normal.");
  fpCommand->SetParameter(parameter);
  parameter  =  new G4UIparameter("ny",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Component of plane normal.");
  fpCommand->SetParameter(parameter);
  parameter  =  new G4UIparameter("nz",'d',omitable = true);
  parameter  -> SetDefaultValue  (0);
  parameter  -> SetGuidance      ("Component of plane normal.");
  fpCommand->SetParameter(parameter);
}

G4VisCommandViewerChangeCutawayPlane::~G4VisCommandViewerChangeCutawayPlane () {
  delete fpCommand;
}

G4String G4VisCommandViewerChangeCutawayPlane::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandViewerChangeCutawayPlane::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
  "ERROR: No current viewer - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }

  size_t index;
  G4double x, y, z, nx, ny, nz;
  G4String unit;
  std::istringstream is (newValue);
  is >> index >> x >> y >> z >> unit >> nx >> ny >> nz;
  G4double F = G4UIcommand::ValueOf(unit);
  x *= F; y *= F; z *= F;

  G4ViewParameters vp = viewer->GetViewParameters();
  vp.ChangeCutawayPlane(index,
			G4Plane3D(G4Normal3D(nx,ny,nz), G4Point3D(x,y,z)));
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Cutaway planes for viewer \"" << viewer->GetName() << "\" now:";
    const G4Planes& cutaways = vp.GetCutawayPlanes();
    for (size_t i = 0; i < cutaways.size(); ++i)
      G4cout << "\n  " << i << ": " << cutaways[i];
    G4cout << G4endl;
  }

  SetViewParameters(viewer, vp);
}

////////////// /vis/viewer/clear ///////////////////////////////////////

G4VisCommandViewerClear::G4VisCommandViewerClear () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/clear", this);
  fpCommand -> SetGuidance ("Clears viewer.");
  fpCommand -> SetGuidance 
    ("By default, clears current viewer.  Specified viewer becomes current."
     "\n\"/vis/viewer/list\" to see  possible viewer names.");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
				 currentAsDefault = true);
}

G4VisCommandViewerClear::~G4VisCommandViewerClear () {
  delete fpCommand;
}

G4String G4VisCommandViewerClear::GetCurrentValue (G4UIcommand*) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  return viewer ? viewer -> GetName () : G4String("none");
}

void G4VisCommandViewerClear::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String& clearName = newValue;
  G4VViewer* viewer = fpVisManager -> GetViewer (clearName);
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Viewer \"" << clearName
	     << "\" not found - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }

  viewer->SetView();
  viewer->ClearView();
  viewer->FinishView();
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Viewer \"" << clearName << "\" cleared." << G4endl;
  }

}

////////////// /vis/viewer/clearCutawayPlanes ///////////////////////////////////////

G4VisCommandViewerClearCutawayPlanes::G4VisCommandViewerClearCutawayPlanes () {
  fpCommand = new G4UIcmdWithoutParameter
    ("/vis/viewer/clearCutawayPlanes", this);
  fpCommand -> SetGuidance ("Clear cutaway planes of current viewer.");
}

G4VisCommandViewerClearCutawayPlanes::~G4VisCommandViewerClearCutawayPlanes () {
  delete fpCommand;
}

G4String G4VisCommandViewerClearCutawayPlanes::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandViewerClearCutawayPlanes::SetNewValue (G4UIcommand*, G4String) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
  "ERROR: No current viewer - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }

  G4ViewParameters vp = viewer->GetViewParameters();
  vp.ClearCutawayPlanes();
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Cutaway planes for viewer \"" << viewer->GetName()
	   << "\" now cleared." << G4endl;
  }

  SetViewParameters(viewer, vp);
}

////////////// /vis/viewer/clearTransients //////////////////////////

G4VisCommandViewerClearTransients::G4VisCommandViewerClearTransients () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/clearTransients", this);
  fpCommand -> SetGuidance ("Clears transients from viewer.");
  fpCommand -> SetGuidance 
  ("By default, operates on current viewer.  Specified viewer becomes current."
     "\n\"/vis/viewer/list\" to see  possible viewer names.");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
				 currentAsDefault = true);
}

G4VisCommandViewerClearTransients::~G4VisCommandViewerClearTransients () {
  delete fpCommand;
}

G4String G4VisCommandViewerClearTransients::GetCurrentValue (G4UIcommand*) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  return viewer ? viewer -> GetName () : G4String("none");
}

void G4VisCommandViewerClearTransients::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String& clearName = newValue;
  G4VViewer* viewer = fpVisManager -> GetViewer (clearName);
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Viewer \"" << clearName
	     << "\" not found - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }

  G4VSceneHandler* sceneHandler = viewer->GetSceneHandler();
  sceneHandler->SetMarkForClearingTransientStore(false);
  fpVisManager->ResetTransientsDrawnFlags();
  sceneHandler->ClearTransientStore();
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Viewer \"" << clearName << "\" cleared of transients."
	   << G4endl;
  }

}

////////////// /vis/viewer/clearVisAttributesModifiers ///////////////////////////////////////

G4VisCommandViewerClearVisAttributesModifiers::G4VisCommandViewerClearVisAttributesModifiers () {
  fpCommand = new G4UIcmdWithoutParameter
  ("/vis/viewer/clearVisAttributesModifiers", this);
  fpCommand -> SetGuidance ("Clear vis attribute modifiers of current viewer.");
  fpCommand -> SetGuidance ("(These are used for touchables, etc.)");
}

G4VisCommandViewerClearVisAttributesModifiers::~G4VisCommandViewerClearVisAttributesModifiers () {
  delete fpCommand;
}

G4String G4VisCommandViewerClearVisAttributesModifiers::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandViewerClearVisAttributesModifiers::SetNewValue (G4UIcommand*, G4String) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
      "ERROR: No current viewer - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }

  G4ViewParameters vp = viewer->GetViewParameters();
  vp.ClearVisAttributesModifiers();
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Vis attributes modifiers for viewer \"" << viewer->GetName()
	   << "\" now cleared." << G4endl;
  }

  SetViewParameters(viewer, vp);
}

////////////// /vis/viewer/clone ///////////////////////////////////////

G4VisCommandViewerClone::G4VisCommandViewerClone () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/viewer/clone", this);
  fpCommand -> SetGuidance ("Clones viewer.");
  fpCommand -> SetGuidance 
    ("By default, clones current viewer.  Clone becomes current."
     "\nClone name, if not provided, is derived from the original name."
     "\n\"/vis/viewer/list\" to see  possible viewer names.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("original-viewer-name", 's', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("clone-name", 's', omitable = true);
  parameter -> SetDefaultValue ("none");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandViewerClone::~G4VisCommandViewerClone () {
  delete fpCommand;
}

G4String G4VisCommandViewerClone::GetCurrentValue (G4UIcommand*) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  G4String originalName = viewer ? viewer -> GetName () : G4String("none");
  return "\"" + originalName + "\"";
}

void G4VisCommandViewerClone::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String originalName, cloneName;
  std::istringstream is (newValue);

  // Need to handle the possibility that the names contain embedded
  // blanks within quotation marks...
  char c = ' ';
  while (is.get(c) && c == ' '){}
  if (c == '"') {
    while (is.get(c) && c != '"') {originalName += c;}
  }
  else {
    originalName += c;
    while (is.get(c) && c != ' ') {originalName += c;}
  }
  originalName = originalName.strip (G4String::both, ' ');
  originalName = originalName.strip (G4String::both, '"');

  G4VViewer* originalViewer = fpVisManager -> GetViewer (originalName);
  if (!originalViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Viewer \"" << originalName
	     << "\" not found - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }
  originalName = originalViewer->GetName();  // Ensures long name.

  while (is.get(c) && c == ' '){}
  if (c == '"') {
    while (is.get(c) && c != '"') {cloneName += c;}
  }
  else {
    cloneName += c;
    while (is.get(c) && c != ' ') {cloneName += c;}
  }
  cloneName = cloneName.strip (G4String::both, ' ');
  cloneName = cloneName.strip (G4String::both, '"');

  G4bool errorWhileNaming = false;
  if (cloneName == "none") {
    G4int subID = 0;
    do {
      cloneName = originalName;
      std::ostringstream oss;
      oss << '-' << subID++;
      G4String::size_type lastDashPosition, nextSpacePosition;
      if ((lastDashPosition = cloneName.rfind('-')) !=  G4String::npos &&
	  (nextSpacePosition = cloneName.find(" ", lastDashPosition)) !=
	  G4String::npos) {
	cloneName.insert(nextSpacePosition, oss.str());
      } else {
	G4String::size_type spacePosition = cloneName.find(' ');
	if (spacePosition != G4String::npos)
	  cloneName.insert(spacePosition, oss.str());
	else
	  errorWhileNaming = true;
      }
    } while (!errorWhileNaming && fpVisManager -> GetViewer (cloneName));
  }

  if (errorWhileNaming) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: While naming clone viewer \"" << cloneName
	     << "\"."
	     << G4endl;
    }
    return;
  }

  if (fpVisManager -> GetViewer (cloneName)) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Putative clone viewer \"" << cloneName
	     << "\" already exists."
	     << G4endl;
    }
    return;
  }

  G4String windowSizeHint =
    originalViewer->GetViewParameters().GetXGeometryString();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  G4int newVerbose(0);
  if (keepVerbose >= 2 ||
      fpVisManager->GetVerbosity() >= G4VisManager::confirmations)
    newVerbose = 2;
  UImanager->SetVerboseLevel(newVerbose);
  UImanager->ApplyCommand(G4String("/vis/viewer/select " + originalName));
  UImanager->ApplyCommand
    (G4String("/vis/viewer/create ! \"" + cloneName + "\" " + windowSizeHint));
  UImanager->ApplyCommand(G4String("/vis/viewer/set/all " + originalName));
  UImanager->SetVerboseLevel(keepVerbose);

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Viewer \"" << originalName << "\" cloned." << G4endl;
    G4cout << "Clone \"" << cloneName << "\" now current." << G4endl;
  }
}

////////////// /vis/viewer/copyViewFrom //////////////////////////

G4VisCommandViewerCopyViewFrom::G4VisCommandViewerCopyViewFrom () {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/copyViewFrom", this);
  fpCommand -> SetGuidance
  ("Copy the camera-specific parameters from the specified viewer.");
  fpCommand -> SetGuidance
  ("Note: To copy ALL view parameters, including scene modifications,"
   "\nuse \"/vis/viewer/set/all\"");
  fpCommand -> SetParameterName ("from-viewer-name", omitable = false);
}

G4VisCommandViewerCopyViewFrom::~G4VisCommandViewerCopyViewFrom () {
  delete fpCommand;
}

G4String G4VisCommandViewerCopyViewFrom::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandViewerCopyViewFrom::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* currentViewer = fpVisManager->GetCurrentViewer();
  if (!currentViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
  "ERROR: G4VisCommandsViewerCopyViewFrom::SetNewValue: no current viewer."
             << G4endl;
    }
    return;
  }

  const G4String& fromViewerName = newValue;
  G4VViewer* fromViewer = fpVisManager -> GetViewer (fromViewerName);
  if (!fromViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Viewer \"" << fromViewerName
	     << "\" not found - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }

  if (fromViewer == currentViewer) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout <<
        "WARNING: G4VisCommandsViewerSet::SetNewValue:"
        "\n  from-viewer and current viewer are identical."
             << G4endl;
    }
    return;
  }

  // Copy camera-specific view parameters
  G4ViewParameters vp = currentViewer->GetViewParameters();
  const G4ViewParameters& fromVP = fromViewer->GetViewParameters();
  vp.SetViewpointDirection  (fromVP.GetViewpointDirection());
  vp.SetLightpointDirection (fromVP.GetLightpointDirection());
  vp.SetLightsMoveWithCamera(fromVP.GetLightsMoveWithCamera());
  vp.SetUpVector            (fromVP.GetUpVector());
  vp.SetFieldHalfAngle      (fromVP.GetFieldHalfAngle());
  vp.SetZoomFactor          (fromVP.GetZoomFactor());
  vp.SetScaleFactor         (fromVP.GetScaleFactor());
  vp.SetCurrentTargetPoint  (fromVP.GetCurrentTargetPoint());
  vp.SetDolly               (fromVP.GetDolly());
  SetViewParameters(currentViewer, vp);
  
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Camera parameters of viewer \"" << currentViewer->GetName()
	   << "\"\n  set to those of viewer \"" << fromViewer->GetName()
	   << "\"."
	   << G4endl;
  }
}

////////////// /vis/viewer/create ///////////////////////////////////////

G4VisCommandViewerCreate::G4VisCommandViewerCreate (): fId (0) {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/viewer/create", this);
  fpCommand -> SetGuidance
    ("Creates a viewer for the specified scene handler.");
  fpCommand -> SetGuidance
    ("Default scene handler is the current scene handler.  Invents a name"
     "\nif not supplied.  (Note: the system adds information to the name"
     "\nfor identification - only the characters up to the first blank are"
     "\nused for removing, selecting, etc.)  This scene handler and viewer"
     "\nbecome current.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("scene-handler", 's', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("viewer-name", 's', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("window-size-hint", 's', omitable = true);
  parameter->SetGuidance
    ("integer (pixels) for square window placed by window manager or"
     " X-Windows-type geometry string, e.g. 600x600-100+100");
  parameter->SetDefaultValue("600");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandViewerCreate::~G4VisCommandViewerCreate () {
  delete fpCommand;
}

G4String G4VisCommandViewerCreate::NextName () {
  std::ostringstream oss;
  G4VSceneHandler* sceneHandler = fpVisManager -> GetCurrentSceneHandler ();
  oss << "viewer-" << fId << " (";
  if (sceneHandler) {
    oss << sceneHandler -> GetGraphicsSystem () -> GetName ();
  }
  else {
    oss << "no_scene_handlers";
  }
  oss << ")";
  return oss.str();
}

G4String G4VisCommandViewerCreate::GetCurrentValue (G4UIcommand*) {
  G4String currentValue;
  G4VSceneHandler* currentSceneHandler =
    fpVisManager -> GetCurrentSceneHandler ();
  if (currentSceneHandler) {
    currentValue = currentSceneHandler -> GetName ();
  }
  else {
    currentValue = "none";
  }
  currentValue += ' ';
  currentValue += '"';
  currentValue += NextName ();
  currentValue += '"';

  currentValue += " 600";  // Default number of pixels for window size hint.

  return currentValue;
}

void G4VisCommandViewerCreate::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String sceneHandlerName, newName;
  G4String windowSizeHintString;
  std::istringstream is (newValue);
  is >> sceneHandlerName;

  // Now need to handle the possibility that the second string
  // contains embedded blanks within quotation marks...
  char c = ' ';
  while (is.get(c) && c == ' '){}
  if (c == '"') {
    while (is.get(c) && c != '"') {newName += c;}
  }
  else {
    newName += c;
    while (is.get(c) && c != ' ') {newName += c;}
  }
  newName = newName.strip (G4String::both, ' ');
  newName = newName.strip (G4String::both, '"');

  // Now get window size hint...
  is >> windowSizeHintString;

  const G4SceneHandlerList& sceneHandlerList =
    fpVisManager -> GetAvailableSceneHandlers ();
  G4int nHandlers = sceneHandlerList.size ();
  if (nHandlers <= 0) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
	"ERROR: G4VisCommandViewerCreate::SetNewValue: no scene handlers."
	"\n  Create a scene handler with \"/vis/sceneHandler/create\""
	     << G4endl;
    }
    return;
  }

  G4int iHandler;
  for (iHandler = 0; iHandler < nHandlers; iHandler++) {
    if (sceneHandlerList [iHandler] -> GetName () == sceneHandlerName) break;
  }

  if (iHandler < 0 || iHandler >= nHandlers) {
    // Invalid command line argument or none.
    // This shouldn't happen!!!!!!
    if (verbosity >= G4VisManager::errors) {
      G4cout << "G4VisCommandViewerCreate::SetNewValue:"
	" invalid scene handler specified."
	    << G4endl;
    }
    return;
  }

  // Valid index.  Set current scene handler and graphics system in
  // preparation for creating viewer.
  G4VSceneHandler* sceneHandler = sceneHandlerList [iHandler];
  if (sceneHandler != fpVisManager -> GetCurrentSceneHandler ()) {
    fpVisManager -> SetCurrentSceneHandler (sceneHandler);
  }

  // Now deal with name of viewer.
  G4String nextName = NextName ();
  if (newName == "") {
    newName = nextName;
  }
  if (newName == nextName) fId++;
  G4String newShortName = fpVisManager -> ViewerShortName (newName);

  for (G4int ih = 0; ih < nHandlers; ih++) {
    G4VSceneHandler* sh = sceneHandlerList [ih];
    const G4ViewerList& viewerList = sh -> GetViewerList ();
    for (size_t iViewer = 0; iViewer < viewerList.size (); iViewer++) {
      if (viewerList [iViewer] -> GetShortName () == newShortName ) {
	if (verbosity >= G4VisManager::errors) {
	  G4cerr << "ERROR: Viewer \"" << newShortName << "\" already exists."
		 << G4endl;
	}
	return;
      }
    }
  }

  // WindowSizeHint and XGeometryString are picked up from the vis
  // manager in the G4VViewer constructor. In G4VisManager, after Viewer
  // creation, we will store theses parameters in G4ViewParameters.

  fpVisManager -> CreateViewer (newName,windowSizeHintString);

  G4VViewer* newViewer = fpVisManager -> GetCurrentViewer ();
  if (newViewer && newViewer -> GetName () == newName) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "New viewer \"" << newName << "\" created." << G4endl;
    }
  }
  else {
    if (verbosity >= G4VisManager::errors) {
      if (newViewer) {
	G4cerr << "ERROR: New viewer doesn\'t match!!!  Curious!!" << G4endl;
      } else {
	G4cout << "WARNING: No viewer created." << G4endl;
      }
    }
  }
  // Refresh if appropriate...
  if (newViewer) {
    if (newViewer->GetViewParameters().IsAutoRefresh()) {
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
    }
    else {
      if (verbosity >= G4VisManager::warnings) {
	G4cout << "Issue /vis/viewer/refresh or flush to see effect."
	       << G4endl;
      }
    }
  }
}

////////////// /vis/viewer/dolly and dollyTo ////////////////////////////

G4VisCommandViewerDolly::G4VisCommandViewerDolly ():
  fDollyIncrement  (0.),
  fDollyTo (0.)
{
  G4bool omitable, currentAsDefault;

  fpCommandDolly = new G4UIcmdWithADoubleAndUnit
    ("/vis/viewer/dolly", this);
  fpCommandDolly -> SetGuidance
    ("Incremental dolly.");
  fpCommandDolly -> SetGuidance
    ("Moves the camera incrementally towards target point.");
  fpCommandDolly -> SetParameterName("increment",
				     omitable=true,
				     currentAsDefault=true);
  fpCommandDolly -> SetDefaultUnit("m");

  fpCommandDollyTo = new G4UIcmdWithADoubleAndUnit
    ("/vis/viewer/dollyTo", this);
  fpCommandDollyTo -> SetGuidance
    ("Dolly to specific coordinate.");
  fpCommandDollyTo -> SetGuidance
 ("Places the camera towards target point relative to standard camera point.");
  fpCommandDollyTo -> SetParameterName("distance",
				       omitable=true,
				       currentAsDefault=true);
  fpCommandDollyTo -> SetDefaultUnit("m");
}

G4VisCommandViewerDolly::~G4VisCommandViewerDolly () {
  delete fpCommandDolly;
  delete fpCommandDollyTo;
}

G4String G4VisCommandViewerDolly::GetCurrentValue (G4UIcommand* command) {
  G4String currentValue;
  if (command == fpCommandDolly) {
    currentValue = fpCommandDolly->ConvertToString(fDollyIncrement, "m");
  }
  else if (command == fpCommandDollyTo) {
    currentValue = fpCommandDollyTo->ConvertToString(fDollyTo, "m");
  }
  return currentValue;
}

void G4VisCommandViewerDolly::SetNewValue (G4UIcommand* command,
					   G4String newValue) {


  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* currentViewer = fpVisManager->GetCurrentViewer();
  if (!currentViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
	"ERROR: G4VisCommandsViewerDolly::SetNewValue: no current viewer."
	     << G4endl;
    }
    return;
  }

  G4ViewParameters vp = currentViewer->GetViewParameters();

  if (command == fpCommandDolly) {
    fDollyIncrement = fpCommandDolly->GetNewDoubleValue(newValue);
    vp.IncrementDolly(fDollyIncrement);
  }
  else if (command == fpCommandDollyTo) {
    fDollyTo = fpCommandDolly->GetNewDoubleValue(newValue);
    vp.SetDolly(fDollyTo);
  }

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Dolly distance changed to " << vp.GetDolly() << G4endl;
  }

  SetViewParameters(currentViewer, vp);
}

////////////// /vis/viewer/flush ///////////////////////////////////////

G4VisCommandViewerFlush::G4VisCommandViewerFlush () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/flush", this);
  fpCommand -> SetGuidance
    ("Compound command: \"/vis/viewer/refresh\" + \"/vis/viewer/update\".");
  fpCommand -> SetGuidance
    ("Useful for refreshing and initiating post-processing for graphics"
     "\nsystems which need post-processing.  By default, acts on current"
     "\nviewer.  \"/vis/viewer/list\" to see possible viewers.  Viewer"
     "\nbecomes current.");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
				 currentAsDefault = true);
}

G4VisCommandViewerFlush::~G4VisCommandViewerFlush () {
  delete fpCommand;
}

G4String G4VisCommandViewerFlush::GetCurrentValue 
(G4UIcommand*) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  return viewer ? viewer -> GetName () : G4String("none");
}

void G4VisCommandViewerFlush::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String& flushName = newValue;
  G4VViewer* viewer = fpVisManager -> GetViewer (flushName);
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Viewer \"" << flushName << "\"" <<
	" not found - \"/vis/viewer/list\"\n  to see possibilities."
	     << G4endl;
    }
    return;
  }

  G4UImanager* ui = G4UImanager::GetUIpointer();
  G4int keepVerbose = ui->GetVerboseLevel();
  G4int newVerbose(0);
  if (keepVerbose >= 2 || verbosity >= G4VisManager::confirmations)
    newVerbose = 2;
  ui->SetVerboseLevel(newVerbose);
  ui->ApplyCommand(G4String("/vis/viewer/refresh " + flushName));
  ui->ApplyCommand(G4String("/vis/viewer/update " + flushName));
  ui->SetVerboseLevel(keepVerbose);
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Viewer \"" << viewer -> GetName () << "\""
	   << " flushed." << G4endl;
  }
}

////////////// /vis/viewer/interpolate ///////////////////////////////////////

G4VisCommandViewerInterpolate::G4VisCommandViewerInterpolate () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/viewer/interpolate", this);
  fpCommand -> SetGuidance
  ("Interpolate views defined by the first argument, which can contain "
   "Unix-shell-style pattern matching characters such as '*', '?' and '[' "
   "- see \"man sh\" and look for \"Pattern Matching\". The contents "
   "of each file are assumed to be \"/vis/viewer\" commands "
   "that specify a particular view. The files are processed in alphanumeric "
   "order of filename. The files may be written by hand or produced by the "
   "\"/vis/viewer/save\" command.");
  fpCommand -> SetGuidance
  ("The default is to search the working directory for files with a .g4view "
   "extension. Another procedure is to assemble view files in a subdirectory, "
   "e.g., \"myviews\"; then they can be interpolated with\n"
   "\"/vis/viewer/interpolate myviews/*\".");
  fpCommand -> SetGuidance
  ("To export interpolated views to file for a future possible movie, "
   "write \"export\" as 5th parameter (OpenGL only).");
  G4UIparameter* parameter;
  parameter = new G4UIparameter("pattern", 's', omitable = true);
  parameter -> SetGuidance("Pattern that defines the view files.");
  parameter -> SetDefaultValue("*.g4view");
  fpCommand -> SetParameter(parameter);
  parameter = new G4UIparameter("no-of-points", 'i', omitable = true);
  parameter -> SetGuidance ("Number of interpolation points per interval.");
  parameter -> SetDefaultValue(50);
  fpCommand -> SetParameter(parameter);
  parameter = new G4UIparameter("wait-time", 's', omitable = true);
  parameter -> SetGuidance("Wait time per interpolated point");
  parameter -> SetDefaultValue("20.");
  fpCommand -> SetParameter(parameter);
  parameter = new G4UIparameter("time-unit", 's', omitable = true);
  parameter -> SetDefaultValue("millisecond");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter("export", 's', omitable = true);
  parameter -> SetDefaultValue("no");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandViewerInterpolate::~G4VisCommandViewerInterpolate () {
  delete fpCommand;
}

G4String G4VisCommandViewerInterpolate::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandViewerInterpolate::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* currentViewer = fpVisManager->GetCurrentViewer();
  if (!currentViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
      "ERROR: G4VisCommandViewerInterpolate::SetNewValue: no current viewer."
	     << G4endl;
    }
    return;
  }

  G4String pattern;
  G4int nInterpolationPoints;
  G4String waitTimePerPointString;
  G4String timeUnit;
  G4String exportString;

  std::istringstream iss (newValue);
  iss
  >> pattern
  >> nInterpolationPoints
  >> waitTimePerPointString
  >> timeUnit
  >> exportString;
  G4String waitTimePerPointDimString(waitTimePerPointString + ' ' + timeUnit);
  const G4double waitTimePerPoint =
  G4UIcommand::ConvertToDimensionedDouble(waitTimePerPointDimString.c_str());
  G4int waitTimePerPointmilliseconds = waitTimePerPoint/millisecond;
  if (waitTimePerPointmilliseconds < 0) waitTimePerPointmilliseconds = 0;

  G4UImanager* uiManager = G4UImanager::GetUIpointer();

  // Save current view parameters
  G4ViewParameters saveVP = currentViewer->GetViewParameters();

  // Save current verbosities
  G4VisManager::Verbosity keepVisVerbosity = fpVisManager->GetVerbosity();
  G4int keepUIVerbosity = uiManager->GetVerboseLevel();

  // Set verbosities for this operation
  fpVisManager->SetVerboseLevel(G4VisManager::errors);
  uiManager->SetVerboseLevel(0);

  // Switch off auto-refresh while we read in the view files (it will be
  // restored later).  Note: the view files do not set auto-refresh.
  G4ViewParameters non_auto = saveVP;
  non_auto.SetAutoRefresh(false);
  currentViewer->SetViewParameters(non_auto);

  // View vector of way points
  std::vector<G4ViewParameters> viewVector;

  const G4int safety = 9999;
  G4int safetyCount = 0;
  G4String pathname;

#ifndef WIN32

  // Execute pattern and get resulting list of filles
  G4String shellCommand = "echo " + pattern;
  FILE *filelist = popen(shellCommand.c_str(), "r");
  if (!filelist) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr
      << "ERROR: G4VisCommandViewerInterpolate::SetNewValue:"
      << "\n  Error obtaining pipe."
	     << G4endl;
    }
    return;
  }

  // Build view vector of way points
  const size_t BUFLENGTH = 999999;
  char buf[BUFLENGTH];
  fgets(buf, BUFLENGTH, filelist);
  std::istringstream fileliststream(buf);
  while (fileliststream >> pathname
         && safetyCount++ < safety) {  // Loop checking, 16.02.2016, J.Allison
    uiManager->ApplyCommand("/control/execute " + pathname);
    viewVector.push_back(currentViewer->GetViewParameters());
  }
  pclose(filelist);

#else // WIN32 (popen is not available in Windows)

  std::experimental::filesystem::v1::path filePattern(pattern);

  // Default pattern : *.g4view
  // Translated to a regexp : ^.*\\.g4view
  // Convert pattern into a regexp
  std::string regexp_pattern("^" + filePattern.filename().string());
  std::string result_pattern = "";
  // Replace '.' by "\\."
  size_t currentPos = 0;
  size_t nextPos = 0;
  std::string currentReplacement = "";
  size_t pos1 = regexp_pattern.find('.', nextPos);
  size_t pos2 = regexp_pattern.find('*', nextPos);
  size_t pos3 = regexp_pattern.find('?', nextPos);
  while ((pos1 != std::string::npos) || (pos2 != std::string::npos) || (pos3 != std::string::npos)) {
    nextPos = pos1;
    currentReplacement = "\\.";
    if (pos2 < nextPos) {
      nextPos = pos2;
      currentReplacement = ".*";
    }
    if (pos3 < nextPos) {
      nextPos = pos3;
      currentReplacement = "(.{1,1})";
    }
    result_pattern += regexp_pattern.substr(currentPos, nextPos - currentPos) + currentReplacement;
    nextPos++;
    currentPos = nextPos;
    pos1 = regexp_pattern.find('.', currentPos);
    pos2 = regexp_pattern.find('*', currentPos);
    pos3 = regexp_pattern.find('?', currentPos);
  }
  result_pattern += regexp_pattern.substr(currentPos);

  // Build view vector of way points
  // Add "./" for empty paths
  G4String parentPath(filePattern.parent_path().string().length() ? filePattern.parent_path().string() : std::string("./"));
  // Iterate through files in directory and apply regex match to filter appropriate files
  std::regex result_pattern_regex (result_pattern, std::regex_constants::basic | std::regex_constants::icase);
  for (auto iter = std::experimental::filesystem::v1::directory_iterator(parentPath);
       iter != std::experimental::filesystem::v1::directory_iterator() && safetyCount++ < safety;
       ++iter)
  {
    const auto& file = iter->path();

    G4String filename(file.filename().string());
    if (std::regex_match(filename, result_pattern_regex))
    {
      uiManager->ApplyCommand("/control/execute " + filename);
      viewVector.push_back(currentViewer->GetViewParameters());
    }
  }

#endif // WIN32

  if (safetyCount >= safety) {
    if (verbosity >= G4VisManager::errors) {
      G4cout <<
      "/vis/viewer/interpolate:"
      "\n  the number of way points exceeds the maximum currently allowed: "
      << safety << G4endl;
    }
    return;
  }

  // Interpolate views
  safetyCount = 0;
  do {
    G4ViewParameters* vp =
    G4ViewParameters::CatmullRomCubicSplineInterpolation(viewVector,nInterpolationPoints);
    if (!vp) break;  // Finished.
    // Set original auto-refresh status.
    vp->SetAutoRefresh(saveVP.IsAutoRefresh());
    currentViewer->SetViewParameters(*vp);
    currentViewer->RefreshView();
    if (exportString == "export" &&
        currentViewer->GetName().contains("OpenGL")) {
      uiManager->ApplyCommand("/vis/ogl/export");
    }
    // File-writing viewers need to close the file
    currentViewer->ShowView();
#ifdef G4VIS_USE_STD11
    if (waitTimePerPointmilliseconds > 0)
      std::this_thread::sleep_for(std::chrono::milliseconds(waitTimePerPointmilliseconds));
#endif
  } while (safetyCount++ < safety);  // Loop checking, 16.02.2016, J.Allison

  // Restore original verbosities
  uiManager->SetVerboseLevel(keepUIVerbosity);
  fpVisManager->SetVerboseLevel(keepVisVerbosity);

  // Restore original view parameters
  currentViewer->SetViewParameters(saveVP);
  currentViewer->RefreshView();
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Viewer \"" << currentViewer -> GetName () << "\""
    << " restored." << G4endl;
  }
}

////////////// /vis/viewer/list ///////////////////////////////////////

G4VisCommandViewerList::G4VisCommandViewerList () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/viewer/list", this);
  fpCommand -> SetGuidance ("Lists viewers(s).");
  fpCommand -> SetGuidance
    ("See \"/vis/verbose\" for definition of verbosity.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter("viewer-name", 's',
				omitable = true);
  parameter -> SetDefaultValue ("all");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("verbosity", 's',
				 omitable = true);
  parameter -> SetDefaultValue ("warnings");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandViewerList::~G4VisCommandViewerList () {
  delete fpCommand;
}

G4String G4VisCommandViewerList::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandViewerList::SetNewValue (G4UIcommand*, G4String newValue) {
  G4String name, verbosityString;
  std::istringstream is (newValue);
  is >> name >> verbosityString;
  G4String shortName = fpVisManager -> ViewerShortName (name);
  G4VisManager::Verbosity verbosity =
    fpVisManager->GetVerbosityValue(verbosityString);

  const G4VViewer* currentViewer = fpVisManager -> GetCurrentViewer ();
  G4String currentViewerShortName;
  if (currentViewer) {
    currentViewerShortName = currentViewer -> GetShortName ();
  }
  else {
    currentViewerShortName = "none";
  }

  const G4SceneHandlerList& sceneHandlerList =
    fpVisManager -> GetAvailableSceneHandlers ();
  G4int nHandlers = sceneHandlerList.size ();
  G4bool found = false;
  G4bool foundCurrent = false;
  for (int iHandler = 0; iHandler < nHandlers; iHandler++) {
    G4VSceneHandler* sceneHandler = sceneHandlerList [iHandler];
    const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
    G4cout
    << "Scene handler \"" << sceneHandler -> GetName () << "\" ("
    << sceneHandler->GetGraphicsSystem()->GetNickname() << ')';
    const G4Scene* pScene = sceneHandler -> GetScene ();
    if (pScene) {
      G4cout << ", scene \"" << pScene -> GetName () << "\"";
    }
    G4cout << ':';
    G4int nViewers = viewerList.size ();
    if (nViewers == 0) {
      G4cout << "\n            No viewers for this scene handler." << G4endl;
    }
    else {
      for (int iViewer = 0; iViewer < nViewers; iViewer++) {
	const G4VViewer* thisViewer = viewerList [iViewer];
	G4String thisName = thisViewer -> GetName ();
	G4String thisShortName = thisViewer -> GetShortName ();
	if (name != "all") {
	  if (thisShortName != shortName) continue;
	}
	found = true;
	G4cout << "\n  ";
	if (thisShortName == currentViewerShortName) {
	  foundCurrent = true;
	  G4cout << "(current)";
	}
	else {
	  G4cout << "         ";
	}
	G4cout << " viewer \"" << thisName << "\"";
	if (verbosity >= G4VisManager::parameters) {
	  G4cout << "\n  " << *thisViewer;
	}
      }
    }
    G4cout << G4endl;
  }

  if (!foundCurrent) {
    G4cout << "No valid current viewer - please create or select one."
	   << G4endl;
  }

  if (!found) {
    G4cout << "No viewers";
    if (name != "all") {
      G4cout << " of name \"" << name << "\"";
    }
    G4cout << " found." << G4endl;
  }
}

////////////// /vis/viewer/pan and panTo ////////////////////////////

G4VisCommandViewerPan::G4VisCommandViewerPan ():
  fPanIncrementRight  (0.),
  fPanIncrementUp  (0.),
  fPanToRight (0.),
  fPanToUp (0.)
{
  G4bool omitable;

  fpCommandPan = new G4UIcommand
    ("/vis/viewer/pan", this);
  fpCommandPan -> SetGuidance
    ("Incremental pan.");
  fpCommandPan -> SetGuidance
    ("Moves the camera incrementally right and up by these amounts (as seen"
     "\nfrom viewpoint direction).");
  G4UIparameter* parameter;
  parameter = new G4UIparameter("right-increment", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandPan -> SetParameter (parameter);
  parameter = new G4UIparameter("up-increment", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandPan -> SetParameter (parameter);
  parameter = new G4UIparameter ("unit", 's', omitable = true);
  parameter -> SetDefaultValue ("m");
  fpCommandPan -> SetParameter (parameter);

  fpCommandPanTo = new G4UIcommand
    ("/vis/viewer/panTo", this);
  fpCommandPanTo -> SetGuidance
    ("Pan to specific coordinate.");
  fpCommandPanTo -> SetGuidance
    ("Places the camera in this position right and up relative to standard"
     "\ntarget point (as seen from viewpoint direction).");
  parameter = new G4UIparameter("right", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandPanTo -> SetParameter (parameter);
  parameter = new G4UIparameter("up", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommandPanTo -> SetParameter (parameter);
  parameter = new G4UIparameter ("unit", 's', omitable = true);
  parameter -> SetDefaultValue ("m");
  fpCommandPanTo -> SetParameter (parameter);
}

G4VisCommandViewerPan::~G4VisCommandViewerPan () {
  delete fpCommandPan;
  delete fpCommandPanTo;
}

G4String G4VisCommandViewerPan::GetCurrentValue (G4UIcommand* command) {
  G4String currentValue;
  if (command == fpCommandPan) {
    currentValue = ConvertToString(fPanIncrementRight, fPanIncrementUp, "m");
  }
  else if (command == fpCommandPanTo) {
    currentValue = ConvertToString(fPanToRight, fPanToUp, "m");
  }
  return currentValue;
}

void G4VisCommandViewerPan::SetNewValue (G4UIcommand* command,
					 G4String newValue) {


  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* currentViewer = fpVisManager->GetCurrentViewer();
  if (!currentViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
	"ERROR: G4VisCommandsViewerPan::SetNewValue: no current viewer."
	     << G4endl;
    }
    return;
  }

  G4ViewParameters vp = currentViewer->GetViewParameters();

  if (command == fpCommandPan) {
    ConvertToDoublePair(newValue, fPanIncrementRight, fPanIncrementUp);
    vp.IncrementPan(fPanIncrementRight, fPanIncrementUp);
  }
  else if (command == fpCommandPanTo) {
    ConvertToDoublePair(newValue, fPanToRight, fPanToUp);
    vp.SetPan(fPanToRight, fPanToUp);
  }

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Current target point now " << vp.GetCurrentTargetPoint()
	   << G4endl;
  }

  SetViewParameters(currentViewer, vp);
}

////////////// /vis/viewer/rebuild ///////////////////////////////////////

G4VisCommandViewerRebuild::G4VisCommandViewerRebuild () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/rebuild", this);
  fpCommand -> SetGuidance ("Forces rebuild of graphical database.");
  fpCommand -> SetGuidance 
    ("By default, acts on current viewer.  \"/vis/viewer/list\""
     "\nto see possible viewers.  Viewer becomes current.");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
				 currentAsDefault = true);
}

G4VisCommandViewerRebuild::~G4VisCommandViewerRebuild () {
  delete fpCommand;
}

G4String G4VisCommandViewerRebuild::GetCurrentValue (G4UIcommand*) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (viewer) {
    return viewer -> GetName ();
  }
  else {
    return "none";
  }
}

void G4VisCommandViewerRebuild::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String& rebuildName = newValue;

  G4VViewer* viewer = fpVisManager -> GetViewer (rebuildName);
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Viewer \"" << rebuildName
	     << "\" not found - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }

  G4VSceneHandler* sceneHandler = viewer->GetSceneHandler();
  if (!sceneHandler) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Viewer \"" << viewer->GetName() << "\"" <<
	" has no scene handler - report serious bug."
	     << G4endl;
    }
    return;
  }

  sceneHandler->ClearTransientStore();
  viewer->NeedKernelVisit();
  viewer->SetView();
  viewer->ClearView();
  viewer->DrawView();

  // Check auto-refresh and print confirmations.
  RefreshIfRequired(viewer);
}

////////////// /vis/viewer/refresh ///////////////////////////////////////

G4VisCommandViewerRefresh::G4VisCommandViewerRefresh () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/refresh", this);
  fpCommand -> SetGuidance
    ("Refreshes viewer.");
  fpCommand -> SetGuidance 
    ("By default, acts on current viewer.  \"/vis/viewer/list\""
     "\nto see possible viewers.  Viewer becomes current.");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
				 currentAsDefault = true);
}

G4VisCommandViewerRefresh::~G4VisCommandViewerRefresh () {
  delete fpCommand;
}

G4String G4VisCommandViewerRefresh::GetCurrentValue (G4UIcommand*) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  return viewer ? viewer -> GetName () : G4String("none");
}

void G4VisCommandViewerRefresh::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4String& refreshName = newValue;
  G4VViewer* viewer = fpVisManager -> GetViewer (refreshName);
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Viewer \"" << refreshName << "\"" <<
	" not found - \"/vis/viewer/list\"\n  to see possibilities."
	     << G4endl;
    }
    return;
  }

  G4VSceneHandler* sceneHandler = viewer->GetSceneHandler();
  if (!sceneHandler) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Viewer \"" << refreshName << "\"" <<
	" has no scene handler - report serious bug."
	     << G4endl;
    }
    return;
  }

  G4Scene* scene = sceneHandler->GetScene();
  if (!scene) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "NOTE: SceneHandler \"" << sceneHandler->GetName()
	     << "\", to which viewer \"" << refreshName << "\"" <<
	"\n  is attached, has no scene - \"/vis/scene/create\" and"
	" \"/vis/sceneHandler/attach\""
	"\n  (or use compound command \"/vis/drawVolume\")."
	     << G4endl;
    }
    return;
  }
  if (scene->GetRunDurationModelList().empty()) {
    G4bool successful = scene -> AddWorldIfEmpty (warn);
    if (!successful) {
      if (verbosity >= G4VisManager::warnings) {
	G4cout <<
	  "WARNING: Scene is empty.  Perhaps no geometry exists."
	  "\n  Try /run/initialize."
	       << G4endl;
      }
      return;
    }
    // Scene has changed.  UpdateVisManagerScene issues
    // /vis/scene/notifyHandlers, which does a refresh anyway, so the
    // ordinary refresh becomes part of the else phrase...
    UpdateVisManagerScene(scene->GetName());
  } else {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Refreshing viewer \"" << viewer -> GetName () << "\"..."
	     << G4endl;
    }
    viewer -> SetView ();
    viewer -> ClearView ();
    viewer -> DrawView ();
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Viewer \"" << viewer -> GetName () << "\"" << " refreshed."
	"\n  (You might also need \"/vis/viewer/update\".)" << G4endl;
    }
  }
}

////////////// /vis/viewer/reset ///////////////////////////////////////

G4VisCommandViewerReset::G4VisCommandViewerReset () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/reset", this);
  fpCommand -> SetGuidance ("Resets viewer.");
  fpCommand -> SetGuidance 
    ("By default, acts on current viewer.  \"/vis/viewer/list\""
     "\nto see possible viewers.  Viewer becomes current.");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
				 currentAsDefault = true);
}

G4VisCommandViewerReset::~G4VisCommandViewerReset () {
  delete fpCommand;
}

G4String G4VisCommandViewerReset::GetCurrentValue (G4UIcommand*) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (viewer) {
    return viewer -> GetName ();
  }
  else {
    return "none";
  }
}

void G4VisCommandViewerReset::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String& resetName = newValue;
  G4VViewer* viewer = fpVisManager -> GetViewer (resetName);
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Viewer \"" << resetName
	     << "\" not found - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }

  viewer->ResetView();
  RefreshIfRequired(viewer);
}

////////////// /vis/viewer/save ///////////////////////////////////////

G4VisCommandViewerSave::G4VisCommandViewerSave () {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/save", this);
  fpCommand -> SetGuidance
  ("Write commands that define the current view to file.");
  fpCommand -> SetGuidance
  ("Read them back into the same or any viewer with \"/control/execute\".");
  fpCommand -> SetGuidance
  ("If the filename is omitted the view is saved to a file "
   "\"g4_nn.g4view\", where nn is a sequential two-digit number.");
  fpCommand -> SetGuidance
  ("If the filename is \"-\", the data are written to G4cout.");
  fpCommand -> SetGuidance
  ("If you are wanting to save views for future interpolation a recommended "
   "procedure is: save views to \"g4_nn.g4view\", as above, then move the files "
   "into a sub-directory, say, \"views\", then interpolate with"
   "\"/vis/viewer/interpolate views/\" (note the trailing \'/\').");
  fpCommand -> SetParameterName ("filename", omitable = true);
  fpCommand -> SetDefaultValue ("");
}

G4VisCommandViewerSave::~G4VisCommandViewerSave () {
  delete fpCommand;
}

G4String G4VisCommandViewerSave::GetCurrentValue
(G4UIcommand*) {
  return "";
}

namespace {
  void WriteCommands
  (std::ostream& os,
   const G4ViewParameters& vp,
   const G4Point3D& stp)  // Standard Target Point
  {
    os
    << vp.CameraAndLightingCommands(stp)
    << vp.DrawingStyleCommands()
    << vp.SceneModifyingCommands()
    << vp.TouchableCommands()
    << std::endl;
  }
}

void G4VisCommandViewerSave::SetNewValue (G4UIcommand*, G4String newValue) {
  
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  
  const G4VViewer* currentViewer = fpVisManager->GetCurrentViewer();
  if (!currentViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
      "ERROR: G4VisCommandsViewerSave::SetNewValue: no current viewer."
      << G4endl;
    }
    return;
  }
  
  const G4Scene* currentScene = currentViewer->GetSceneHandler()->GetScene();
  if (!currentScene) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
      "ERROR: G4VisCommandsViewerSave::SetNewValue: no current scene."
      << G4endl;
    }
    return;
  }

  // Get view parameters and ther relevant information.
  G4ViewParameters vp = currentViewer->GetViewParameters();
  // Concatenate any private vis attributes modifiers...
  const std::vector<G4ModelingParameters::VisAttributesModifier>*
  privateVAMs = currentViewer->GetPrivateVisAttributesModifiers();
  if (privateVAMs) {
    std::vector<G4ModelingParameters::VisAttributesModifier>::const_iterator i;
    for (i = privateVAMs->begin(); i != privateVAMs->end(); ++i) {
      vp.AddVisAttributesModifier(*i);
    }
  }
  const G4Point3D& stp = currentScene->GetStandardTargetPoint();

  G4String filename = newValue;

  if (newValue.length() == 0) {
    // Null filename - generate a filename
    const G4int maxNoOfFiles = 100;
    static G4int sequenceNumber = 0;
    if (sequenceNumber >= maxNoOfFiles) {
      if (verbosity >= G4VisManager::errors) {
        G4cerr
        << "ERROR: G4VisCommandsViewerSave::SetNewValue: Maximum number, "
        << maxNoOfFiles
        << ", of files exceeded."
        << G4endl;
      }
      return;
    }
    std::ostringstream oss;
    oss << std::setw(2) << std::setfill('0') << sequenceNumber++;
    filename = "g4_" + oss.str() + ".g4view";
  }

  if (filename == "-") {
    // Write to standard output
    WriteCommands(G4cout,vp,stp);
  } else {
    // Write to file - but add extension if not prescribed
    if (!filename.contains('.')) {
      // No extension supplied - add .g4view
      filename += ".g4view";
    }
    std::ofstream ofs(filename);
    if (!ofs) {
      if (verbosity >= G4VisManager::errors) {
        G4cerr <<
        "ERROR: G4VisCommandsViewerSave::SetNewValue: Trouble opening file \""
        << filename << "\"."
        << G4endl;
      }
      ofs.close();
      return;
    }
    WriteCommands(ofs,vp,stp);
    ofs.close();
  }
  
  if (verbosity >= G4VisManager::warnings) {
    G4cout << "Viewer \"" << currentViewer -> GetName ()
    << "\"" << " saved to ";
    if (filename == "-") {
      G4cout << "G4cout.";
    } else {
      G4cout << "file \'" << filename << "\"." <<
      "\n  Read the view back into this or any viewer with"
      "\n  \"/control/execute " << filename << "\" or use"
      "\n  \"/vis/viewer/interpolate\" if you have several saved files -"
      "\n  see \"help /vis/viewer/interpolate\" for guidance.";
    }
    G4cout << G4endl;
  }
}

////////////// /vis/viewer/scale and scaleTo ////////////////////////////

G4VisCommandViewerScale::G4VisCommandViewerScale ():
  fScaleMultiplier (G4Vector3D (1., 1., 1.)),
  fScaleTo         (G4Vector3D (1., 1., 1.))
{
  G4bool omitable, currentAsDefault;

  fpCommandScale = new G4UIcmdWith3Vector
    ("/vis/viewer/scale", this);
  fpCommandScale -> SetGuidance ("Incremental (non-uniform) scaling.");
  fpCommandScale -> SetGuidance
    ("Multiplies components of current scaling by components of this factor."
     "\n Scales (x,y,z) by corresponding components of the resulting factor.");
  fpCommandScale -> SetGuidance
    ("");
  fpCommandScale -> SetParameterName
    ("x-scale-multiplier","y-scale-multiplier","z-scale-multiplier",
     omitable=true, currentAsDefault=true);

  fpCommandScaleTo = new G4UIcmdWith3Vector
    ("/vis/viewer/scaleTo", this);
  fpCommandScaleTo -> SetGuidance ("Absolute (non-uniform) scaling.");
  fpCommandScaleTo -> SetGuidance
    ("Scales (x,y,z) by corresponding components of this factor.");
  fpCommandScaleTo -> SetParameterName
    ("x-scale-factor","y-scale-factor","z-scale-factor",
     omitable=true, currentAsDefault=true);
}

G4VisCommandViewerScale::~G4VisCommandViewerScale () {
  delete fpCommandScale;
  delete fpCommandScaleTo;
}

G4String G4VisCommandViewerScale::GetCurrentValue (G4UIcommand* command) {
  G4String currentValue;
  if (command == fpCommandScale) {
    currentValue = fpCommandScale->ConvertToString(G4ThreeVector(fScaleMultiplier));
  }
  else if (command == fpCommandScaleTo) {
    currentValue = fpCommandScaleTo->ConvertToString(G4ThreeVector(fScaleTo));
  }
  return currentValue;
}

void G4VisCommandViewerScale::SetNewValue (G4UIcommand* command,
					   G4String newValue) {


  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* currentViewer = fpVisManager->GetCurrentViewer();
  if (!currentViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
	"ERROR: G4VisCommandsViewerScale::SetNewValue: no current viewer."
	     << G4endl;
    }
    return;
  }

  G4ViewParameters vp = currentViewer->GetViewParameters();

  if (command == fpCommandScale) {
    fScaleMultiplier = fpCommandScale->GetNew3VectorValue(newValue);
    vp.MultiplyScaleFactor(fScaleMultiplier);
  }
  else if (command == fpCommandScaleTo) {
    fScaleTo = fpCommandScale->GetNew3VectorValue(newValue);
    vp.SetScaleFactor(fScaleTo);
  }

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Scale factor changed to " << vp.GetScaleFactor() << G4endl;
  }

  SetViewParameters(currentViewer, vp);
}

////////////// /vis/viewer/select ///////////////////////////////////////

G4VisCommandViewerSelect::G4VisCommandViewerSelect () {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/select", this);
  fpCommand -> SetGuidance ("Selects viewer.");
  fpCommand -> SetGuidance
    ("Specify viewer by name.  \"/vis/viewer/list\" to see possible viewers.");
  fpCommand -> SetParameterName ("viewer-name", omitable = false);
}

G4VisCommandViewerSelect::~G4VisCommandViewerSelect () {
  delete fpCommand;
}

G4String G4VisCommandViewerSelect::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandViewerSelect::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String& selectName = newValue;
  G4VViewer* viewer = fpVisManager -> GetViewer (selectName);

  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Viewer \"" << selectName << "\"";
      G4cerr << " not found - \"/vis/viewer/list\""
	"\n  to see possibilities."
	     << G4endl;
    }
    return;
  }

  if (viewer == fpVisManager -> GetCurrentViewer ()) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout << "WARNING: Viewer \"" << viewer -> GetName () << "\""
	     << " already selected." << G4endl;
    }
    return;
  }

  // Set pointers, call SetView and print confirmation.
  fpVisManager -> SetCurrentViewer (viewer);

  RefreshIfRequired(viewer);
}

////////////// /vis/viewer/update ///////////////////////////////////////

G4VisCommandViewerUpdate::G4VisCommandViewerUpdate () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/update", this);
  fpCommand -> SetGuidance
    ("Triggers graphical database post-processing for viewers"
     "\nusing that technique.");
  fpCommand -> SetGuidance
    ("For such viewers the view only becomes visible with this command."
     "\nBy default, acts on current viewer.  \"/vis/viewer/list\""
     "\nto see possible viewers.  Viewer becomes current.");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
				 currentAsDefault = true);
}

G4VisCommandViewerUpdate::~G4VisCommandViewerUpdate () {
  delete fpCommand;
}

G4String G4VisCommandViewerUpdate::GetCurrentValue (G4UIcommand*) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (viewer) {
    return viewer -> GetName ();
  }
  else {
    return "none";
  }
}

void G4VisCommandViewerUpdate::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String& updateName = newValue;

  G4VViewer* viewer = fpVisManager -> GetViewer (updateName);
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cout <<
	"WARNING: command \"/vis/viewer/update\" could not be applied: no current viewer."
	     << G4endl;
    }
    return;
  }

  G4VSceneHandler* sceneHandler = viewer->GetSceneHandler();
  if (!sceneHandler) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Viewer \"" << updateName << "\"" <<
	" has no scene handler - report serious bug."
	     << G4endl;
    }
    return;
  }

  G4Scene* scene = sceneHandler->GetScene();
  if (!scene) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "NOTE: SceneHandler \"" << sceneHandler->GetName()
	     << "\", to which viewer \"" << updateName << "\"" <<
	"\n  is attached, has no scene - \"/vis/scene/create\" and"
	" \"/vis/sceneHandler/attach\""
	"\n  (or use compound command \"/vis/drawVolume\")."
	     << G4endl;
    }
    return;
  }

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Viewer \"" << viewer -> GetName () << "\"";
    G4cout << " post-processing triggered." << G4endl;
  }
  viewer -> ShowView ();
  // Assume future need to "refresh" transients...
  sceneHandler -> SetMarkForClearingTransientStore(true);
}

////////////// /vis/viewer/zoom and zoomTo ////////////////////////////

G4VisCommandViewerZoom::G4VisCommandViewerZoom ():
  fZoomMultiplier (1.),
  fZoomTo         (1.)
{
  G4bool omitable, currentAsDefault;

  fpCommandZoom = new G4UIcmdWithADouble
    ("/vis/viewer/zoom", this);
  fpCommandZoom -> SetGuidance ("Incremental zoom.");
  fpCommandZoom -> SetGuidance
    ("Multiplies current magnification by this factor.");
  fpCommandZoom -> SetParameterName("multiplier",
				     omitable=true,
				     currentAsDefault=true);

  fpCommandZoomTo = new G4UIcmdWithADouble
    ("/vis/viewer/zoomTo", this);
  fpCommandZoomTo -> SetGuidance ("Absolute zoom.");
  fpCommandZoomTo -> SetGuidance
    ("Magnifies standard magnification by this factor.");
  fpCommandZoomTo -> SetParameterName("factor",
				       omitable=true,
				       currentAsDefault=true);
}

G4VisCommandViewerZoom::~G4VisCommandViewerZoom () {
  delete fpCommandZoom;
  delete fpCommandZoomTo;
}

G4String G4VisCommandViewerZoom::GetCurrentValue (G4UIcommand* command) {
  G4String currentValue;
  if (command == fpCommandZoom) {
    currentValue = fpCommandZoom->ConvertToString(fZoomMultiplier);
  }
  else if (command == fpCommandZoomTo) {
    currentValue = fpCommandZoomTo->ConvertToString(fZoomTo);
  }
  return currentValue;
}

void G4VisCommandViewerZoom::SetNewValue (G4UIcommand* command,
					  G4String newValue) {


  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* currentViewer = fpVisManager->GetCurrentViewer();
  if (!currentViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
	"ERROR: G4VisCommandsViewerZoom::SetNewValue: no current viewer."
	     << G4endl;
    }
    return;
  }

  G4ViewParameters vp = currentViewer->GetViewParameters();

  if (command == fpCommandZoom) {
    fZoomMultiplier = fpCommandZoom->GetNewDoubleValue(newValue);
    vp.MultiplyZoomFactor(fZoomMultiplier);
  }
  else if (command == fpCommandZoomTo) {
    fZoomTo = fpCommandZoom->GetNewDoubleValue(newValue);
    vp.SetZoomFactor(fZoomTo);
  }

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Zoom factor changed to " << vp.GetZoomFactor() << G4endl;
  }

  SetViewParameters(currentViewer, vp);
}
