// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsViewer.cc,v 1.10 2000-01-11 17:22:31 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/viewer commands - John Allison  25th October 1998

#include "G4VisCommandsViewer.hh"

#include "G4VisManager.hh"
#include "G4GraphicsSystemList.hh"
#include "G4VisCommandsScene.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"
#include "g4std/strstream"

G4VVisCommandViewer::G4VVisCommandViewer () {}

G4VVisCommandViewer::~G4VVisCommandViewer () {}

G4String G4VVisCommandViewer::ShortName (const G4String& name) {
  G4String shortName (name);
  if (shortName.contains (' ')) {
    shortName = shortName (0, shortName.first (' '));
  }
  return shortName;
}

void G4VVisCommandViewer::UpdateCandidateLists () {

  const G4SceneHandlerList& sceneHandlerList =
    fpVisManager -> GetAvailableSceneHandlers ();
  G4int nHandlers = sceneHandlerList.entries ();

  G4String viewerNameList;
  for (int iHandler = 0; iHandler < nHandlers; iHandler++) {
    G4VSceneHandler* sceneHandler = sceneHandlerList [iHandler];
    const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
    for (int iViewer = 0; iViewer < viewerList.entries (); iViewer++) {
      const G4String& viewerName = viewerList [iViewer] -> GetName ();
      viewerNameList += ShortName (viewerName) + " ";
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

  const G4SceneHandlerList& sceneHandlerList = fpVisManager -> GetAvailableSceneHandlers ();
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

  for (iHandler = 0; iHandler < nHandlers; iHandler++) {
    G4VSceneHandler* sceneHandler = sceneHandlerList [iHandler];
    const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
    for (int iViewer = 0; iViewer < viewerList.entries (); iViewer++) {
      if (viewerList [iViewer] -> GetName () == newName ) {
	G4cout << "Viewer \"" << newName << "\" already exists." << G4endl;
	return;
      }
    }
  }

  // Create viewer.
  fpVisManager -> CreateViewer (newName);
  G4cout << "New viewer \"" << newName << "\" created." << G4endl;

  UpdateCandidateLists ();
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
  G4String shortName = ShortName (name);

  const G4VViewer* currentViewer = fpVisManager -> GetCurrentViewer ();
  G4String currentViewerShortName;
  if (currentViewer) {
    currentViewerShortName = ShortName (currentViewer -> GetName ());
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
	G4String thisShortName = ShortName (thisName);
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
    G4cout << "No valid current viewer - please create or select one." << G4endl;
  }

  if (!found) {
    G4cout << "No viewers";
    if (name != "all") {
      G4cout << " of name \"" << name << "\"";
    }
    G4cout << " found." << endl;
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
  G4String refreshShortName = ShortName (refreshName);

  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  G4bool found = true;
  if (refreshShortName != ShortName (viewer -> GetName ())) {
    found = false;
    const G4SceneHandlerList& sceneHandlerList =
      fpVisManager -> GetAvailableSceneHandlers ();
    G4int nHandlers = sceneHandlerList.entries ();
    for (G4int iHandler = 0; iHandler < nHandlers; iHandler++) {
      G4VSceneHandler* sceneHandler = sceneHandlerList [iHandler];
      const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
      for (G4int iViewer = 0; iViewer < viewerList.entries (); iViewer++) {
	viewer = viewerList [iViewer];
	if (refreshShortName == ShortName (viewer -> GetName ())) {
	  found = true;
	  fpVisManager -> SetCurrentViewer (viewer);
	  break;
	}
      }
      if (found) break;
    }
  }

  G4cout << "Viewer \"" << refreshName << "\"";
  if (found) {
    viewer -> ClearView ();
    viewer -> SetViewParameters (fpVisManager -> GetCurrentViewParameters ());
    viewer -> DrawView ();
    viewer -> ShowView ();
    G4cout << " refreshed.";
  }
  else {
    G4cout << " not found - \"/vis/viewer/list\""
      "\n  to see possibilities.";
  }
  G4cout << endl;
}

////////////// /vis/viewer/remove ///////////////////////////////////////

G4VisCommandViewerRemove::G4VisCommandViewerRemove () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/remove", this);
  fpCommand -> SetGuidance ("/vis/viewer/remove <viewer-name>");
  fpCommand -> SetGuidance ("Removes viewers.");
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
  G4String removeShortName = ShortName (removeName);
  G4VViewer* currentViewer = fpVisManager -> GetCurrentViewer ();
  G4String currentName;
  if (currentViewer) {
    currentName = currentViewer -> GetName ();
  }
  else {
    currentName = "none";
  }
  G4String currentShortName = ShortName (currentName);

  const G4SceneHandlerList& sceneHandlerList = fpVisManager -> GetAvailableSceneHandlers ();
  G4int nHandlers = sceneHandlerList.entries ();
  G4int iHandler, iViewer;
  G4VSceneHandler* sceneHandler;
  G4VViewer* viewer;
  G4bool found = false;
  for (iHandler = 0; iHandler < nHandlers; iHandler++) {
    sceneHandler = sceneHandlerList [iHandler];
    const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
    for (iViewer = 0; iViewer < viewerList.entries (); iViewer++) {
      viewer = viewerList [iViewer];
      const G4String& viewerName = viewer -> GetName ();
      if (removeShortName == ShortName (viewerName)) {
	found = true;
	break;
      }
    }
  }

  if (!found) {
    G4cout << "Viewer \"" << removeName
	   << "\" not found - \"/vis/viewer/list\" to see possibilities."
	   << G4endl;
      return;
  }

  G4cout << "Viewer \"" << removeName << "\" removed." << G4endl;
  if (removeShortName == currentShortName) {
    fpVisManager -> DeleteCurrentViewer ();
  }
  else {
    sceneHandler -> SetViewerList ().remove (viewer);
    G4cout << "Current viewer unchanged." << G4endl;
  }

  UpdateCandidateLists ();
}

////////////// /vis/viewer/select ///////////////////////////////////////

G4VisCommandViewerSelect::G4VisCommandViewerSelect () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/select", this);
  fpCommand -> SetGuidance ("/vis/viewer/select [<viewer-name>]");
  fpCommand -> SetGuidance ("Selects current viewer.");
  fpCommand -> SetGuidance
    ("Specify viewer by name (\"/vis/viewer/list\""
     "\n  to see possibilities).");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
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
  G4String selectShortName = ShortName (selectName);

  const G4SceneHandlerList& sceneHandlerList =
    fpVisManager -> GetAvailableSceneHandlers ();
  G4int nHandlers = sceneHandlerList.entries ();
  G4int iHandler, iViewer;
  G4VSceneHandler* sceneHandler;
  G4VViewer* viewer;
  G4bool found = false;
  for (iHandler = 0; iHandler < nHandlers; iHandler++) {
    sceneHandler = sceneHandlerList [iHandler];
    const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
    for (iViewer = 0; iViewer < viewerList.entries (); iViewer++) {
      viewer = viewerList [iViewer];
      if (selectShortName == ShortName (viewer -> GetName ())) {
	found = true;
	break;
      }
    }
    if (found) break;
  }

  G4cout << "Viewer \"" << selectName << "\"";
  if (found) {
    G4VViewer* currentViewer = fpVisManager -> GetCurrentViewer ();
    if (currentViewer &&
	ShortName (currentViewer -> GetName ())	== selectShortName) {
      G4cout << " already selected." << G4endl;
    }
    else {
      G4cout << " being selected." << G4endl;
      fpVisManager -> SetCurrentViewer (viewer);
    }
  }
  else {
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
  G4String showShortName = ShortName (showName);

  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  G4bool found = true;
  if (showShortName != ShortName (viewer -> GetName ())) {
    found = false;
    const G4SceneHandlerList& sceneHandlerList =
      fpVisManager -> GetAvailableSceneHandlers ();
    G4int nHandlers = sceneHandlerList.entries ();
    for (G4int iHandler = 0; iHandler < nHandlers; iHandler++) {
      G4VSceneHandler* sceneHandler = sceneHandlerList [iHandler];
      const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
      for (G4int iViewer = 0; iViewer < viewerList.entries (); iViewer++) {
	viewer = viewerList [iViewer];
	if (showShortName == ShortName (viewer -> GetName ())) {
	  found = true;
	  fpVisManager -> SetCurrentViewer (viewer);
	  break;
	}
      }
      if (found) break;
    }
  }

  G4cout << "Viewer \"" << showName << "\"";
  if (found) {
    viewer -> ShowView ();
    G4cout << " post-processing triggered.";
  }
  else {
    G4cout << " not found - \"/vis/viewer/list\""
      "\n  to see possibilities.";
  }
  G4cout << G4endl;
}
