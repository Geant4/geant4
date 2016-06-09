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
// $Id: G4VisCommandsViewer.cc,v 1.53 2005/11/22 17:00:15 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $

// /vis/viewer commands - John Allison  25th October 1998

#include "G4VisCommandsViewer.hh"

#include "G4VisManager.hh"
#include "G4GraphicsSystemList.hh"
#include "G4VisCommandsScene.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include <sstream>

G4VVisCommandViewer::G4VVisCommandViewer () {}

G4VVisCommandViewer::~G4VVisCommandViewer () {}

void G4VVisCommandViewer::SetViewParameters
(G4VViewer* viewer, const G4ViewParameters& viewParams) {
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  viewer->SetViewParameters(viewParams);
  G4VSceneHandler* sceneHandler = viewer->GetSceneHandler();
  if (sceneHandler && sceneHandler->GetScene()) {
    if (viewParams.IsAutoRefresh()) {
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
    }
    else {
      if (verbosity >= G4VisManager::confirmations) {
	G4cout << "Issue /vis/viewer/refresh to see effect." << G4endl;
      }
    }
  }
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
      G4cout << "ERROR: Viewer \"" << clearName
	     << "\" not found - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }

  viewer->ClearView();
  viewer->FinishView();
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Viewer \"" << clearName << "\" cleared." << G4endl;
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
  char c;
  while (is.get(c) && c == ' ');
  if (c == '"') {
    while (is.get(c) && c != '"') newName += c;
  }
  else {
    newName += c;
    while (is.get(c) && c != ' ') newName += c;
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
      G4cout <<
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

  for (iHandler = 0; iHandler < nHandlers; iHandler++) {
    G4VSceneHandler* sceneHandler = sceneHandlerList [iHandler];
    const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
    for (size_t iViewer = 0; iViewer < viewerList.size (); iViewer++) {
      if (viewerList [iViewer] -> GetShortName () == newShortName ) {
	if (verbosity >= G4VisManager::errors) {
	  G4cout << "ERROR: Viewer \"" << newShortName << "\" already exists."
		 << G4endl;
	}
	return;
      }
    }
  }

  // Parse windowSizeHintString to extract first field for backwards
  // compatibility...
  std::istringstream issw;
  G4int windowSizeHint;
  size_t i;
  for (i = 0; i < windowSizeHintString.size(); ++i) {
    char c = windowSizeHintString[i];
    if (c == 'x' || c == 'X' || c == '+' || c == '-') break;
  }
  if (i != windowSizeHintString.size()) {
    // x or X or + or - found - must be a X-Window-type geometry string...
    // Pick out the first field for backwards compatibility...
    issw.str(windowSizeHintString.substr(0,i));
    issw >> windowSizeHint;
  } else { // ...old-style integer...
    issw.str(windowSizeHintString);
    if (!(issw >> windowSizeHint)) {
      if (verbosity >= G4VisManager::errors) {
	G4cout << "ERROR: Unrecognised geometry string \""
	       << windowSizeHintString
	       << "\".  Using 600."
	       << G4endl;
      }
      windowSizeHint = 600;
    }
    // Reconstitute windowSizeHintString...
    std::ostringstream ossw;
    ossw << windowSizeHint << 'x' << windowSizeHint;
    windowSizeHintString = ossw.str();
  }
  fpVisManager->SetWindowSizeHint (windowSizeHint, windowSizeHint);
  fpVisManager->SetXGeometryString(windowSizeHintString);
  // WindowSizeHint and XGeometryString are picked up from the vis
  // manager in the G4VViewer constructor.  They have to be held by
  // the vis manager until the viewer is contructed - next line...

  // Create viewer.
  fpVisManager -> CreateViewer (newName);
  G4VViewer* newViewer = fpVisManager -> GetCurrentViewer ();
  if (newViewer && newViewer -> GetName () == newName) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "New viewer \"" << newName << "\" created." << G4endl;
    }
  }
  else {
    if (verbosity >= G4VisManager::errors) {
      if (newViewer) {
	G4cout << "ERROR: New viewer doesn\'t match!!!  Curious!!" << G4endl;
      } else {
	G4cout << "WARNING: No viewer created." << G4endl;
      }
    }
  }
  // Refresh if appropriate...
  if (newViewer) SetViewParameters(newViewer, newViewer->GetViewParameters());
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
      G4cout <<
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
      G4cout << "ERROR: Viewer \"" << flushName << "\"" <<
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
  parameter -> SetDefaultValue (0);
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

  const G4SceneHandlerList& sceneHandlerList = fpVisManager -> GetAvailableSceneHandlers ();
  G4int nHandlers = sceneHandlerList.size ();
  G4bool found = false;
  G4bool foundCurrent = false;
  for (int iHandler = 0; iHandler < nHandlers; iHandler++) {
    G4VSceneHandler* sceneHandler = sceneHandlerList [iHandler];
    const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
    G4cout << "Scene handler \"" << sceneHandler -> GetName ();
    const G4Scene* pScene = sceneHandler -> GetScene ();
    if (pScene) {
      G4cout << "\", scene \"" << pScene -> GetName () << "\":";
    }
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
      G4cout <<
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
      G4cout << "ERROR: Viewer \"" << rebuildName
	     << "\" not found - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }

  viewer->SetNeedKernelVisit(true);
  SetViewParameters(viewer, viewer->GetDefaultViewParameters());
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
      G4cout << "ERROR: Viewer \"" << refreshName << "\"" <<
	" not found - \"/vis/viewer/list\"\n  to see possibilities."
	     << G4endl;
    }
    return;
  }

  G4VSceneHandler* sceneHandler = viewer->GetSceneHandler();
  if (!sceneHandler) {
    if (verbosity >= G4VisManager::errors) {
      G4cout << "ERROR: Viewer \"" << refreshName << "\"" <<
	" has no scene handler - report serious bug."
	     << G4endl;
    }
    return;
  }

  G4Scene* scene = sceneHandler->GetScene();
  if (!scene) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout << "WARNING: SceneHandler \"" << sceneHandler->GetName()
	     << "\", to which viewer \"" << refreshName << "\"" <<
	"\n  is attached, has no scene - \"/vis/scene/create\" and"
	"\"/vis/sceneHandler/attach\""
	"\n  (or use compound command \"/vis/drawVolume\")."
	     << G4endl;
    }
    return;
  }
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
      G4cout << "ERROR: Viewer \"" << resetName
	     << "\" not found - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }

  SetViewParameters(viewer, viewer->GetDefaultViewParameters());
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
      G4cout <<
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
      G4cout << "ERROR: Viewer \"" << selectName << "\"";
      G4cout << " not found - \"/vis/viewer/list\""
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

  fpVisManager -> SetCurrentViewer (viewer);  // Prints confirmation.

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
	"ERROR: G4VisCommandsViewerUpdate::SetNewValue: no current viewer."
	     << G4endl;
    }
    return;
  }

  G4VSceneHandler* sceneHandler = viewer->GetSceneHandler();
  if (!sceneHandler) {
    if (verbosity >= G4VisManager::errors) {
      G4cout << "ERROR: Viewer \"" << updateName << "\"" <<
	" has no scene handler - report serious bug."
	     << G4endl;
    }
    return;
  }

  G4Scene* scene = sceneHandler->GetScene();
  if (!scene) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout << "WARNING: SceneHandler \"" << sceneHandler->GetName()
	     << "\", to which viewer \"" << updateName << "\"" <<
	"\n  is attached, has no scene - \"/vis/scene/create\" and"
	"\"/vis/sceneHandler/attach\""
	"\n  (or use compound command \"/vis/drawVolume\")."
	     << G4endl;
    }
    return;
  }

  if (viewer) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Viewer \"" << viewer -> GetName () << "\"";
      G4cout << " post-processing triggered." << G4endl;
    }
    viewer -> ShowView ();
  }
  else {
    if (verbosity >= G4VisManager::errors) {
      G4cout << "ERROR: Viewer \"" << updateName << "\"";
      G4cout << " not found - \"/vis/viewer/list\""
	"\n  to see possibilities." << G4endl;
    }
  }
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
      G4cout <<
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
