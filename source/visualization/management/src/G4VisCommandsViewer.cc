// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsViewer.cc,v 1.7 1999-11-05 16:31:39 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/viewer commands - John Allison  25th October 1998

#include "G4VisCommandsViewer.hh"

#include "G4VisManager.hh"
#include "G4GraphicsSystemList.hh"
#include "G4VisCommandsScene.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"
#ifdef WIN32
#include <strstrea.h>
#else
#include <strstream.h>
#endif

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

  const G4SceneHandlerList& sceneHandlerList = fpVisManager -> GetAvailableSceneHandlers ();
  G4int nHandlers = sceneHandlerList.entries ();

  fViewerNameList = G4String ();
  for (int iHandler = 0; iHandler < nHandlers; iHandler++) {
    G4VSceneHandler* sceneHandler = sceneHandlerList [iHandler];
    const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
    for (int iViewer = 0; iViewer < viewerList.entries (); iViewer++) {
      const G4String& viewerName = viewerList [iViewer] -> GetName ();
      fViewerNameList += ShortName (viewerName) + " ";
    }
  }
  fViewerNameList = fViewerNameList.strip ();

  if (fpCommandViewerRemove) {
    fpCommandViewerRemove -> GetParameter (0) ->
      SetParameterCandidates (fViewerNameList);
  }

  if (fpCommandViewerSelect) {
    fpCommandViewerSelect -> GetParameter (0) ->
      SetParameterCandidates (fViewerNameList);
  }

  if (fpCommandViewerUpdate) {
    fpCommandViewerUpdate -> GetParameter (0) ->
      SetParameterCandidates (fViewerNameList);
  }
}

G4String G4VVisCommandViewer::fViewerNameList;

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
  fpCommandViewerCreate = fpCommand;
}

G4VisCommandViewerCreate::~G4VisCommandViewerCreate () {
  delete fpCommand;
}

G4String G4VisCommandViewerCreate::NextName () {
  const int charLength = 100;
  char nextName [charLength];
  ostrstream ost (nextName, charLength);
  G4VSceneHandler* sceneHandler = fpVisManager -> GetCurrentSceneHandler ();
  ost << "viewer-" << fId << " (";
  if (sceneHandler) {
    ost << sceneHandler -> GetGraphicsSystem () -> GetName ();
  }
  else {
    ost << "no_scene_handlers";
  }
  ost << ")" << ends;
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
  istrstream is ((char*)newValue.data());
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
	   << endl;
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
	   << endl;
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
	G4cout << "Viewer \"" << newName << "\" already exists." << endl;
	return;
      }
    }
  }

  // Create viewer.
  fpVisManager -> CreateViewer (newName);
  G4cout << "New viewer \"" << newName << "\" created." << endl;

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
  istrstream is ((char*)newValue.data());
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
      G4cout << "\n            No viewers for this scene handler." << endl;
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
      G4cout << endl;
    }
  }

  if (!foundCurrent) {
    G4cout << "No valid current viewer - please create or select one." << endl;
  }

  if (!found) {
    G4cout << "No viewers";
    if (name != "all") {
      G4cout << " of name \"" << name << "\"";
    }
    G4cout << " found." << endl;
  }
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
  fpCommandViewerRemove = fpCommand;
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
	   << endl;
      return;
  }

  G4cout << "Viewer \"" << removeName << "\" removed." << endl;
  if (removeShortName == currentShortName) {
    fpVisManager -> DeleteCurrentViewer ();
  }
  else {
    sceneHandler -> SetViewerList ().remove (viewer);
    G4cout << "Current viewer unchanged." << endl;
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
  fpCommandViewerSelect = fpCommand;
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
      G4cout << " already selected." << endl;
    }
    else {
      G4cout << " being selected." << endl;
      fpVisManager -> SetCurrentViewer (viewer);
    }
  }
  else {
    G4cout << " not found - \"/vis/viewer/list\""
      "\n  to see possibilities."
	   << endl;
  }
}

////////////// /vis/viewer/update ///////////////////////////////////////

G4VisCommandViewerUpdate::G4VisCommandViewerUpdate () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/viewer/update", this);
  fpCommand -> SetGuidance ("/vis/viewer/update [<viewer-name>]");
  fpCommand -> SetGuidance
    ("Updates viewer (necessary for graphical database post-processing.");
  fpCommand -> SetGuidance ("Viewer becomes current.");
  fpCommand -> SetGuidance
    ("Specify viewer by name (\"/vis/viewer/list\""
     "\n  to see possibilities).");
  fpCommand -> SetParameterName ("viewer-name",
				 omitable = true,
				 currentAsDefault = true);
  fpCommandViewerUpdate = fpCommand;
}

G4VisCommandViewerUpdate::~G4VisCommandViewerUpdate () {
  delete fpCommand;
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
  G4String& updateName = newValue;
  G4String updateShortName = ShortName (updateName);

  G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
  G4bool found = true;
  if (updateShortName != ShortName (viewer -> GetName ())) {
    found = false;
    const G4SceneHandlerList& sceneHandlerList =
      fpVisManager -> GetAvailableSceneHandlers ();
    G4int nHandlers = sceneHandlerList.entries ();
    for (G4int iHandler = 0; iHandler < nHandlers; iHandler++) {
      G4VSceneHandler* sceneHandler = sceneHandlerList [iHandler];
      const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
      for (G4int iViewer = 0; iViewer < viewerList.entries (); iViewer++) {
	viewer = viewerList [iViewer];
	if (updateShortName == ShortName (viewer -> GetName ())) {
	  found = true;
	  fpVisManager -> SetCurrentViewer (viewer);
	  fpVisManager -> SetCurrentSceneHandler (sceneHandler);
	  fpVisManager -> SetCurrentGraphicsSystem
	    (sceneHandler -> GetGraphicsSystem ());
	  break;
	}
      }
      if (found) break;
    }
  }

  G4cout << "Viewer \"" << updateName << "\"";
  if (found) {
    viewer -> ShowView ();
    G4cout << " updated.";
  }
  else {
    G4cout << " not found - \"/vis/viewer/list\""
      "\n  to see possibilities.";
  }
  G4cout << endl;
}
