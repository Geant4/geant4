// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsViewer.cc,v 1.16 2000-05-04 19:17:14 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/viewer commands - John Allison  25th October 1998

#include "G4VisCommandsViewer.hh"

#include "G4VisManager.hh"
#include "G4GraphicsSystemList.hh"
#include "G4VisCommandsScene.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "g4std/strstream"

G4VVisCommandViewer::G4VVisCommandViewer () {}

G4VVisCommandViewer::~G4VVisCommandViewer () {}

void G4VVisCommandViewer::UpdateCandidateLists () {

  const G4SceneHandlerList& sceneHandlerList =
    fpVisManager -> GetAvailableSceneHandlers ();
  G4int nHandlers = sceneHandlerList.entries ();

  G4String viewerNameList;
  for (int iHandler = 0; iHandler < nHandlers; iHandler++) {
    G4VSceneHandler* sceneHandler = sceneHandlerList [iHandler];
    const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
    for (int iViewer = 0; iViewer < viewerList.entries (); iViewer++) {
      viewerNameList += viewerList [iViewer] -> GetShortName () + " ";
    }
  }
  viewerNameList = viewerNameList.strip ();
  viewerNameCommandsIterator i;
  for (i = viewerNameCommands.begin (); i != viewerNameCommands.end (); ++i) {
    (*i)->GetParameter (0) -> SetParameterCandidates (viewerNameList);
  }
}

////////////// /vis/viewer/create ///////////////////////////////////////

G4VisCommandViewerCreate::G4VisCommandViewerCreate (): fId (0) {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/viewer/create", this);
  fpCommand -> SetGuidance
    ("/vis/viewer/create  [<scene-handler>] [<viewer-name>]");
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
  fpCommand -> SetGuidance ("This scene handler and viewer become current.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("scene-handler", 's', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("viewer-name", 's', omitable = true);
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

  /* Replaces currentValue += NextName ();
  // Now need to handle the possibility that NextName contains
  // embedded blanks.
  currentValue += '"';
  const int charLength = 100;
  char nextName [charLength];
  strncpy (nextName, NextName (), charLength);
  for (int i = 0; nextName[i]; i++) {
    currentValue += nextName[i];
  }
  */

  return currentValue;
}

void G4VisCommandViewerCreate::SetNewValue (G4UIcommand* command,
					   G4String newValue) {
  G4String sceneHandlerName, newName;
  G4std::istrstream is ((char*)newValue.data());
  is >> sceneHandlerName;
  // Now need to handle the possibility that the second string
  // contains embedded blanks within quotation marks.
  char c;
  while (is.get (c)) newName += c;
  newName = newName.strip (G4String::both, ' ');
  newName = newName.strip (G4String::both, '"');

  const G4SceneHandlerList& sceneHandlerList =
    fpVisManager -> GetAvailableSceneHandlers ();
  G4int nHandlers = sceneHandlerList.entries ();
  if (nHandlers <= 0) {
    G4cout << "G4VisCommandViewerCreate::SetNewValue: no scene handlers."
      "\n  Create a scene handler with \"/vis/sceneHandler/create\""
	   << G4endl;
    return;
  }

  G4int iHandler;
  for (iHandler = 0; iHandler < nHandlers; iHandler++) {
    if (sceneHandlerList [iHandler] -> GetName () == sceneHandlerName) break;
  }

  if (iHandler < 0 || iHandler >= nHandlers) {
    // Invalid command line argument or non.
    // This shouldn't happen!!!!!!
    G4cerr << "G4VisCommandViewerCreate::SetNewValue:"
      " invalid scene handler specified."
	   << G4endl;
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
    for (int iViewer = 0; iViewer < viewerList.entries (); iViewer++) {
      if (viewerList [iViewer] -> GetShortName () == newShortName ) {
	G4cout << "Viewer \"" << newShortName << "\" already exists."
	       << G4endl;
	return;
      }
    }
  }

  // Create viewer.
  fpVisManager -> CreateViewer (newName);
  if (fpVisManager -> GetCurrentViewer () -> GetName () == newName) {
    G4cout << "New viewer \"" << newName << "\" created." << G4endl;
    UpdateCandidateLists ();
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
    ("<verbosity> is 0 for short (default) or 1 for long listing.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter("viewer-name", 's',
				omitable = true);
  parameter -> SetCurrentAsDefault (false);
  parameter -> SetDefaultValue ("all");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("verbosity", 'i',
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
  G4String name;
  G4int verbosity;
  G4std::istrstream is ((char*)newValue.data());
  is >> name >> verbosity;
  G4String shortName = fpVisManager -> ViewerShortName (name);

  const G4VViewer* currentViewer = fpVisManager -> GetCurrentViewer ();
  G4String currentViewerShortName;
  if (currentViewer) {
    currentViewerShortName = currentViewer -> GetShortName ();
  }
  else {
    currentViewerShortName = "none";
  }

  const G4SceneHandlerList& sceneHandlerList = fpVisManager -> GetAvailableSceneHandlers ();
  G4int nHandlers = sceneHandlerList.entries ();
  G4bool found = false;
  G4bool foundCurrent = false;
  for (int iHandler = 0; iHandler < nHandlers; iHandler++) {
    G4VSceneHandler* sceneHandler = sceneHandlerList [iHandler];
    const G4Scene* pScene = sceneHandler -> GetScene ();
    const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
    G4cout << "Scene handler \"" << sceneHandler -> GetName ()
	   << "\", scene \"" << pScene -> GetName () << "\":";
    G4int nViewers = viewerList.entries ();
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
	if (verbosity > 0) {
	  G4cout << "\n  " << *thisViewer << '\n';
	}
      }
      G4cout << G4endl;
    }
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
  G4String& refreshName = newValue;
  G4VViewer* viewer = fpVisManager -> GetViewer (refreshName);
  if (viewer) {
    G4cout << "Refreshing viewer \"" << viewer -> GetName () << "\"..."
	   << G4endl;
    viewer -> ClearView ();
    viewer -> DrawView ();
    G4cout << "Viewer \"" << viewer -> GetName () << "\"" << " refreshed."
	   "\n  (You might also need \"/vis/viewer/show\".)" << G4endl;
  }
  else {
    G4cout << "Viewer \"" << refreshName << "\"" <<
      " not found - \"/vis/viewer/list\"\n  to see possibilities."
	   << G4endl;
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
    G4cout << "Viewer \"" << removeName
	   << "\" not found - \"/vis/viewer/list\" to see possibilities."
	   << G4endl;
      return;
  }

  G4cout << "Viewer \"" << viewer -> GetName () << "\" removed." << G4endl;
  if (viewer -> GetShortName () == currentShortName) {
    fpVisManager -> DeleteCurrentViewer ();
  }
  else {
    G4VSceneHandler* sceneHandler = viewer -> GetSceneHandler ();
    sceneHandler -> SetViewerList ().remove (viewer);
    G4cout << "Current viewer is unchanged (\""
	   << currentViewer -> GetName () << "\")." << G4endl;
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
  G4String& resetName = newValue;
  G4VViewer* viewer = fpVisManager -> GetViewer (resetName);
  if (!viewer) {
    G4cout << "Viewer \"" << resetName
	   << "\" not found - \"/vis/viewer/list\" to see possibilities."
	   << G4endl;
    return;
  }

  if (fpVisManager -> IsValidView ()) {
    G4ViewParameters vp = viewer->GetViewParameters();
    const G4VSceneHandler* sceneHandler = viewer->GetSceneHandler();
    const G4Scene* scene = sceneHandler->GetScene();
    vp.SetCurrentTargetPoint (scene->GetStandardTargetPoint ());
    vp.SetZoomFactor (1.);
    vp.SetDolly (0.);
    vp.SetViewpointDirection (G4Vector3D (0., 0., 1.));
    vp.SetUpVector (G4Vector3D (0., 1., 0.));
    viewer->SetViewParameters(vp);
    G4cout << "Target point reset to centre of scene, "
           << G4BestUnit (scene->GetStandardTargetPoint (), "Length");
    G4cout << ".\nZoom factor reset to 1.";
    G4cout << "\nDolly distance reset to 0.";
    G4cout << "\nViewpoint direction reset to +z.";
    G4cout << "\nUp vector reset to +y.";
    G4cout << "\nViewer \"" << viewer -> GetName () << "\" reset.";
    G4cout << "\nIssue \"/vis/viewer/refresh\" to see effect." << G4endl;
    // For now...
    fpVisManager->SetCurrentViewParameters() = vp;
  }
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
  G4String& selectName = newValue;
  G4VViewer* viewer = fpVisManager -> GetViewer (selectName);

  if (viewer) {
    G4cout << "Viewer \"" << viewer -> GetName () << "\"";
    if (viewer == fpVisManager -> GetCurrentViewer ()) {
      G4cout << " already selected." << G4endl;
    }
    else {
      G4cout << " being selected." << G4endl;
      fpVisManager -> SetCurrentViewer (viewer);
    }
  }
  else {
    G4cout << "Viewer \"" << selectName << "\"";
    G4cout << " not found - \"/vis/viewer/list\""
      "\n  to see possibilities."
	   << G4endl;
  }
}

////////////// /vis/viewer/show ///////////////////////////////////////
// Synonym (deprecated): /vis/viewer/update ///////////////////////////

G4VisCommandViewerShow::G4VisCommandViewerShow () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/show", this);
  fpCommand -> SetGuidance ("/vis/viewer/show [<viewer-name>]");
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
  fpCommand1 = new G4UIcmdWithAString ("/vis/viewer/update", this);
  fpCommand1 -> SetGuidance
    ("Synonym for \"/vis/viewer/show\" - see that command for guidance.");
  fpCommand1 -> SetParameterName ("viewer-name",
				  omitable = true,
				  currentAsDefault = true);
  viewerNameCommands.push_back (fpCommand1);
}

G4VisCommandViewerShow::~G4VisCommandViewerShow () {
  delete fpCommand;
  delete fpCommand1;
}

G4String G4VisCommandViewerShow::GetCurrentValue 
(G4UIcommand* command) {
  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  if (viewer) {
    return viewer -> GetName ();
  }
  else {
    return "none";
  }
}

void G4VisCommandViewerShow::SetNewValue (G4UIcommand* command,
					   G4String newValue) {
  G4String& showName = newValue;
  G4VViewer* viewer = fpVisManager -> GetViewer (showName);

  if (viewer) {
    G4cout << "Viewer \"" << viewer -> GetName () << "\"";
    G4cout << " post-processing triggered." << G4endl;
    viewer -> ShowView ();
  }
  else {
    G4cout << "Viewer \"" << showName << "\"";
    G4cout << " not found - \"/vis/viewer/list\""
      "\n  to see possibilities." << G4endl;
  }
}
