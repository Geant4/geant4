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
// $Id: G4VisCommandsViewer.cc,v 1.35 2002-04-26 21:23:31 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "g4std/strstream"

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

void G4VVisCommandViewer::UpdateCandidateLists () {

  const G4SceneHandlerList& sceneHandlerList =
    fpVisManager -> GetAvailableSceneHandlers ();
  G4int nHandlers = sceneHandlerList.size ();

  G4String viewerNameList;
  for (int iHandler = 0; iHandler < nHandlers; iHandler++) {
    G4VSceneHandler* sceneHandler = sceneHandlerList [iHandler];
    const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
    for (size_t iViewer = 0; iViewer < viewerList.size (); iViewer++) {
      viewerNameList += viewerList [iViewer] -> GetShortName () + " ";
    }
  }
  viewerNameList = viewerNameList.strip ();
  viewerNameCommandsIterator i;
  for (i = viewerNameCommands.begin (); i != viewerNameCommands.end (); ++i) {
    (*i)->GetParameter (0) -> SetParameterCandidates (viewerNameList);
  }
}

////////////// /vis/viewer/clear ///////////////////////////////////////

G4VisCommandViewerClear::G4VisCommandViewerClear () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/clear", this);
  fpCommand -> SetGuidance ("/vis/viewer/clear [<viewer-name>]");
  fpCommand -> SetGuidance ("Clears viewer.");
  fpCommand -> SetGuidance ("Viewer becomes current.");
  fpCommand -> SetGuidance
    ("Specify viewer by name (\"/vis/viewer/list\""
     "\n  to see possibilities).");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
				 currentAsDefault = true);
  viewerNameCommands.push_back (fpCommand);
}

G4VisCommandViewerClear::~G4VisCommandViewerClear () {
  delete fpCommand;
}

G4String G4VisCommandViewerClear::GetCurrentValue (G4UIcommand* command) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (viewer) {
    return viewer -> GetName ();
  }
  else {
    return "none";
  }
}

void G4VisCommandViewerClear::SetNewValue (G4UIcommand* command,
					    G4String newValue) {

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
    ("/vis/viewer/create  [<scene-handler>] [<viewer-name>] [<pixels>]");
  fpCommand -> SetGuidance
    ("Creates an viewer for a specific scene handler.");
  fpCommand -> SetGuidance
    ("Default scene handler is the current scene handler.");
  fpCommand -> SetGuidance
    ("Invents a name if not supplied.  (Note: the system adds information");
  fpCommand -> SetGuidance
    ("to the name for identification - only the characters up to the first");
  fpCommand -> SetGuidance
    ("blank are used for removing, selecting, etc.)");
  fpCommand -> SetGuidance
    ("The 3rd parameter is a window size hint.");
  fpCommand -> SetGuidance ("This scene handler and viewer become current.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("scene-handler", 's', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("viewer-name", 's', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("pixels", 'd', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommand -> SetParameter (parameter);
  viewerNameCommands.push_back (fpCommand);
}

G4VisCommandViewerCreate::~G4VisCommandViewerCreate () {
  delete fpCommand;
}

G4String G4VisCommandViewerCreate::NextName () {
  const int charLength = 100;
  char nextName [charLength];
  G4std::ostrstream ost (nextName, charLength);
  G4VSceneHandler* sceneHandler = fpVisManager -> GetCurrentSceneHandler ();
  ost << "viewer-" << fId << " (";
  if (sceneHandler) {
    ost << sceneHandler -> GetGraphicsSystem () -> GetName ();
  }
  else {
    ost << "no_scene_handlers";
  }
  ost << ")" << G4std::ends;
  return nextName;
}

G4String G4VisCommandViewerCreate::GetCurrentValue (G4UIcommand* command) {
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

void G4VisCommandViewerCreate::SetNewValue (G4UIcommand* command,
					   G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String sceneHandlerName, newName;
  G4int windowSizeHint;
  G4std::istrstream is ((char*)newValue.data());
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

  // Now get number of pixels...
  is >> windowSizeHint;

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
    // Invalid command line argument or non.
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

  fpVisManager->SetWindowSizeHint (windowSizeHint, windowSizeHint);
  // These are picked up in the G4VViewer constructor.  The problem is
  // these have to be set *before* construction, i.e., before we have
  // a viewer.

  // Create viewer.
  fpVisManager -> CreateViewer (newName);
  G4VViewer* newViewer = fpVisManager -> GetCurrentViewer ();
  if (newViewer -> GetName () == newName) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "New viewer \"" << newName << "\" created." << G4endl;
    }
    UpdateCandidateLists ();
  }
  else {
    if (verbosity >= G4VisManager::errors) {
      G4cout << "ERROR: New viewer doesn\'t match!!!  Curious!!" << G4endl;
    }
  }
  // Refresh if appropriate...
  SetViewParameters(newViewer, newViewer->GetViewParameters());
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
    ("/vis/viewer/dolly [<increment>] [<unit>]");
  fpCommandDolly -> SetGuidance
    ("Moves the camera incrementally in by this distance.");
  fpCommandDolly -> SetParameterName("increment",
				     omitable=true,
				     currentAsDefault=true);
  fpCommandDolly -> SetDefaultUnit("m");

  fpCommandDollyTo = new G4UIcmdWithADoubleAndUnit
    ("/vis/viewer/dollyTo", this);
  fpCommandDollyTo -> SetGuidance
    ("/vis/viewer/dollyTo [<distance>] [<unit>]");
  fpCommandDollyTo -> SetGuidance
    ("Moves the camera in this distance relative to standard target point.");
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
  fpCommand -> SetGuidance ("/vis/viewer/flush [<viewer-name>]");
  fpCommand -> SetGuidance
    ("Compound command: /vis/viewer/refresh + /vis/viewer/update.");
  fpCommand -> SetGuidance
    ("Useful for refreshing and initiating post-processing for graphics"
     "\n  systems which need post-processing.  Viewer becomes current.");
  fpCommand -> SetGuidance
    ("Specify viewer by name (\"/vis/viewer/list\""
     "\n  to see possibilities).");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
				 currentAsDefault = true);
  viewerNameCommands.push_back (fpCommand);
}

G4VisCommandViewerFlush::~G4VisCommandViewerFlush () {
  delete fpCommand;
}

G4String G4VisCommandViewerFlush::GetCurrentValue 
(G4UIcommand* command) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (viewer) {
    return viewer -> GetName ();
  }
  else {
    return "none";
  }
}

void G4VisCommandViewerFlush::SetNewValue (G4UIcommand* command,
					   G4String newValue) {

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

//////// /vis/viewer/lightsThetaPhi and lightsVector /////////////

G4VisCommandViewerLights::G4VisCommandViewerLights ():
  fLightsVector (G4ThreeVector(0.,0.,1.))
{
  G4bool omitable;

  fpCommandLightsThetaPhi = new G4UIcommand
    ("/vis/viewer/lightsThetaPhi", this);
  fpCommandLightsThetaPhi -> SetGuidance
    ("/vis/viewer/lightsThetaPhi  [<theta>] [<phi>] [deg|rad]");
  fpCommandLightsThetaPhi -> SetGuidance
    ("Set direction from target to lights.");
  G4UIparameter* parameter;
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
    ("/vis/viewer/lightsVector", this);
  fpCommandLightsVector -> SetGuidance
    ("/vis/viewer/lightsVector  [<x>] [<y>] [<z>]");
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
}

G4VisCommandViewerLights::~G4VisCommandViewerLights () {
  delete fpCommandLightsThetaPhi;
  delete fpCommandLightsVector;
}

G4String G4VisCommandViewerLights::GetCurrentValue (G4UIcommand* command) {
  G4String currentValue;
  if (command == fpCommandLightsThetaPhi) {
    currentValue = ConvertToString(fLightsVector.theta(),
			   fLightsVector.phi(), "deg");
  }
  else if (command == fpCommandLightsVector) {
    currentValue = ConvertToString(fLightsVector);
  }
  return currentValue;
}

void G4VisCommandViewerLights::SetNewValue (G4UIcommand* command,
					   G4String newValue) {
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  if (command == fpCommandLightsThetaPhi) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout <<
	"WARNING: DEPRECATED: use \"/vis/viewer/set/LightsThetaPhi\"."
	     << G4endl;
    }
    G4UImanager::GetUIpointer()->ApplyCommand(
      G4String("/vis/viewer/set/LightsThetaPhi " + newValue));
  }
  else if (command == fpCommandLightsVector) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout <<
	"WARNING: DEPRECATED: use \"/vis/viewer/set/LightsVector\"."
	     << G4endl;
    }
    G4UImanager::GetUIpointer()->ApplyCommand(
      G4String("/vis/viewer/set/LightsVector " + newValue));
  }
}

////////////// /vis/viewer/list ///////////////////////////////////////

G4VisCommandViewerList::G4VisCommandViewerList () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/viewer/list", this);
  fpCommand -> SetGuidance
    ("/vis/viewer/list [<viewer-name>] [<verbosity>]");
  fpCommand -> SetGuidance ("Lists viewers(s).");
  fpCommand -> SetGuidance ("<viewer-name> default is \"all\"");
  fpCommand -> SetGuidance
    ("See /vis/verbose for definition of verbosity.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter("viewer-name", 's',
				omitable = true);
  parameter -> SetCurrentAsDefault (false);
  parameter -> SetDefaultValue ("all");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("verbosity", 's',
				 omitable = true);
  parameter -> SetCurrentAsDefault (false);
  parameter -> SetDefaultValue (0);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandViewerList::~G4VisCommandViewerList () {
  delete fpCommand;
}

G4String G4VisCommandViewerList::GetCurrentValue (G4UIcommand* command) {
  return "";
}

void G4VisCommandViewerList::SetNewValue (G4UIcommand* command,
					 G4String newValue) {
  G4String name, verbosityString;
  G4std::istrstream is ((char*)newValue.data());
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
    ("/vis/viewer/pan [<right-increment>] [<up-increment>] [<unit>]");
  fpCommandPan -> SetGuidance
    ("Moves the camera incrementally right and up by these amounts.");
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
    ("/vis/viewer/panTo [<right>] [<up>] [<unit>]");
  fpCommandPanTo -> SetGuidance
    ("Moves the camera to this position right and up relative to standard"
     "target point (as seen from viewpoint direction).");
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
    GetNewDoublePairValue(newValue, fPanIncrementRight, fPanIncrementUp);
    vp.IncrementPan(fPanIncrementRight, fPanIncrementUp);
  }
  else if (command == fpCommandPanTo) {
    GetNewDoublePairValue(newValue, fPanToRight, fPanToUp);
    vp.SetPan(fPanToRight, fPanToUp);
  }

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Current target point now " << vp.GetCurrentTargetPoint()
	   << G4endl;
  }

  SetViewParameters(currentViewer, vp);
}

////////////// /vis/viewer/refresh ///////////////////////////////////////

G4VisCommandViewerRefresh::G4VisCommandViewerRefresh () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/refresh", this);
  fpCommand -> SetGuidance ("/vis/viewer/refresh [<viewer-name>]");
  fpCommand -> SetGuidance
    ("Refreshes viewer.");
  fpCommand -> SetGuidance ("Viewer becomes current.");
  fpCommand -> SetGuidance
    ("Specify viewer by name (\"/vis/viewer/list\""
     "\n  to see possibilities).");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
				 currentAsDefault = true);
  viewerNameCommands.push_back (fpCommand);
}

G4VisCommandViewerRefresh::~G4VisCommandViewerRefresh () {
  delete fpCommand;
}

G4String G4VisCommandViewerRefresh::GetCurrentValue 
(G4UIcommand* command) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (viewer) {
    return viewer -> GetName ();
  }
  else {
    return "none";
  }
}

void G4VisCommandViewerRefresh::SetNewValue (G4UIcommand* command,
					   G4String newValue) {

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

////////////// /vis/viewer/remove ///////////////////////////////////////

G4VisCommandViewerRemove::G4VisCommandViewerRemove () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/remove", this);
  fpCommand -> SetGuidance ("/vis/viewer/remove <viewer-name>");
  fpCommand -> SetGuidance ("Removes viewer.");
  fpCommand -> SetGuidance
    ("Specify viewer by name (\"/vis/viewer/list\""
     "\n  to see possibilities).");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = false,
				 currentAsDefault = true);
  viewerNameCommands.push_back (fpCommand);
}

G4VisCommandViewerRemove::~G4VisCommandViewerRemove () {
  delete fpCommand;
}

G4String G4VisCommandViewerRemove::GetCurrentValue (G4UIcommand* command) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (viewer) {
    return viewer -> GetName ();
  }
  else {
    return "none";
  }
}

void G4VisCommandViewerRemove::SetNewValue (G4UIcommand* command,
					    G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String& removeName = newValue;

  G4VViewer* currentViewer = fpVisManager -> GetCurrentViewer ();
  G4String currentShortName;
  if (currentViewer) {
    currentShortName = currentViewer -> GetShortName ();
  }
  else {
    currentShortName = "none";
  }

  G4VViewer* viewer = fpVisManager -> GetViewer (removeName);
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cout << "ERROR: Viewer \"" << removeName
	     << "\" not found - \"/vis/viewer/list\" to see possibilities."
	     << G4endl;
    }
    return;
  }

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Viewer \"" << viewer -> GetName () << "\" removed." << G4endl;
  }
  if (viewer -> GetShortName () == currentShortName) {
    fpVisManager -> DeleteCurrentViewer ();
  }
  else {
    G4VSceneHandler* sceneHandler = viewer -> GetSceneHandler ();
    G4ViewerList& viewerList = sceneHandler -> SetViewerList ();
    viewerList.remove (viewer);
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Current viewer is unchanged (\""
	     << currentViewer -> GetName () << "\")." << G4endl;
    }
  }

  UpdateCandidateLists ();
}

////////////// /vis/viewer/reset ///////////////////////////////////////

G4VisCommandViewerReset::G4VisCommandViewerReset () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/reset", this);
  fpCommand -> SetGuidance ("/vis/viewer/reset [<viewer-name>]");
  fpCommand -> SetGuidance ("Resets viewer.");
  fpCommand -> SetGuidance ("Viewer becomes current.");
  fpCommand -> SetGuidance
    ("Specify viewer by name (\"/vis/viewer/list\""
     "\n  to see possibilities).");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
				 currentAsDefault = true);
  viewerNameCommands.push_back (fpCommand);
}

G4VisCommandViewerReset::~G4VisCommandViewerReset () {
  delete fpCommand;
}

G4String G4VisCommandViewerReset::GetCurrentValue (G4UIcommand* command) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (viewer) {
    return viewer -> GetName ();
  }
  else {
    return "none";
  }
}

void G4VisCommandViewerReset::SetNewValue (G4UIcommand* command,
					    G4String newValue) {

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

////////////// /vis/viewer/select ///////////////////////////////////////

G4VisCommandViewerSelect::G4VisCommandViewerSelect () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/select", this);
  fpCommand -> SetGuidance ("/vis/viewer/select <viewer-name>");
  fpCommand -> SetGuidance ("Selects current viewer.");
  fpCommand -> SetGuidance
    ("Specify viewer by name (\"/vis/viewer/list\""
     "\n  to see possibilities).");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = false,
				 currentAsDefault = true);
  viewerNameCommands.push_back (fpCommand);
}

G4VisCommandViewerSelect::~G4VisCommandViewerSelect () {
  delete fpCommand;
}

G4String G4VisCommandViewerSelect::GetCurrentValue 
(G4UIcommand* command) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (viewer) {
    return viewer -> GetName ();
  }
  else {
    return "none";
  }
}

void G4VisCommandViewerSelect::SetNewValue (G4UIcommand* command,
					   G4String newValue) {

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
// Synonym (deprecated): /vis/viewer/show ///////////////////////////

G4VisCommandViewerUpdate::G4VisCommandViewerUpdate () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/update", this);
  fpCommand -> SetGuidance ("/vis/viewer/update [<viewer-name>]");
  fpCommand -> SetGuidance
    ("Triggers graphical database post-processing for viewers"
     " using that technique.");
  fpCommand -> SetGuidance
    ("For such viewers the view only becomes visible with this command.");
  fpCommand -> SetGuidance ("Viewer becomes current.");
  fpCommand -> SetGuidance
    ("Specify viewer by name (\"/vis/viewer/list\""
     "\n  to see possibilities).");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
				 currentAsDefault = true);
  viewerNameCommands.push_back (fpCommand);
  fpCommand1 = new G4UIcmdWithAString ("/vis/viewer/show", this);
  fpCommand1 -> SetGuidance
    ("Synonym for \"/vis/viewer/update\" - see that command for guidance.");
  fpCommand1 -> SetParameterName ("viewer-name",
				  omitable = true,
				  currentAsDefault = true);
  viewerNameCommands.push_back (fpCommand1);
}

G4VisCommandViewerUpdate::~G4VisCommandViewerUpdate () {
  delete fpCommand;
  delete fpCommand1;
}

G4String G4VisCommandViewerUpdate::GetCurrentValue 
(G4UIcommand* command) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (viewer) {
    return viewer -> GetName ();
  }
  else {
    return "none";
  }
}

void G4VisCommandViewerUpdate::SetNewValue (G4UIcommand* command,
					   G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String& updateName = newValue;
  G4VViewer* viewer = fpVisManager -> GetViewer (updateName);

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

//////// /vis/viewer/viewpointThetaPhi and viewpointVector /////////////

G4VisCommandViewerViewpoint::G4VisCommandViewerViewpoint ():
  fViewpointVector (G4ThreeVector(0.,0.,1.))
{
  G4bool omitable;

  fpCommandViewpointThetaPhi = new G4UIcommand
    ("/vis/viewer/viewpointThetaPhi", this);
  fpCommandViewpointThetaPhi -> SetGuidance
    ("/vis/viewer/viewpointThetaPhi  [<theta>] [<phi>] [deg|rad]");
  fpCommandViewpointThetaPhi -> SetGuidance
    ("Set direction from target to camera.  Also changes lightpoint direction"
     "\nif lights are set to move with camera.");
  G4UIparameter* parameter;
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
    ("/vis/viewer/viewpointVector", this);
  fpCommandViewpointVector -> SetGuidance
    ("/vis/viewer/viewpointVector  [<x>] [<y>] [<z>]");
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

G4VisCommandViewerViewpoint::~G4VisCommandViewerViewpoint () {
  delete fpCommandViewpointThetaPhi;
  delete fpCommandViewpointVector;
}

G4String G4VisCommandViewerViewpoint::GetCurrentValue (G4UIcommand* command) {
  G4String currentValue;
  if (command == fpCommandViewpointThetaPhi) {
    currentValue = ConvertToString(fViewpointVector.theta(),
			   fViewpointVector.phi(), "deg");
  }
  else if (command == fpCommandViewpointVector) {
    currentValue = ConvertToString(fViewpointVector);
  }
  return currentValue;
}

void G4VisCommandViewerViewpoint::SetNewValue (G4UIcommand* command,
					   G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  if (command == fpCommandViewpointThetaPhi) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout <<
	"WARNING: DEPRECATED: use \"/vis/viewer/set/viewpointThetaPhi\"."
	     << G4endl;
    }
    G4UImanager::GetUIpointer()->ApplyCommand(
      G4String("/vis/viewer/set/viewpointThetaPhi " + newValue));
  }
  else if (command == fpCommandViewpointVector) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout <<
	"WARNING: DEPRECATED: use \"/vis/viewer/set/viewpointVector\"."
	     << G4endl;
    }
    G4UImanager::GetUIpointer()->ApplyCommand(
      G4String("/vis/viewer/set/viewpointVector " + newValue));
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
  fpCommandZoom -> SetGuidance
    ("/vis/viewer/zoom [<multiplier>]");
  fpCommandZoom -> SetGuidance
    ("Multiplies magnification by this factor.");
  fpCommandZoom -> SetParameterName("multiplier",
				     omitable=true,
				     currentAsDefault=true);

  fpCommandZoomTo = new G4UIcmdWithADouble
    ("/vis/viewer/zoomTo", this);
  fpCommandZoomTo -> SetGuidance
    ("/vis/viewer/zoomTo [<factor>]");
  fpCommandZoomTo -> SetGuidance
    ("Magnifies by this factor relative to standard view.");
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
