// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsScene.cc,v 1.2 1999-01-09 16:31:20 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/scene commands - John Allison  9th August 1998

#include "G4VisCommandsScene.hh"

#include "G4VisManager.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4ApplicationState.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"
#ifdef WIN32
#include <strstrea.h>
#else
#include <strstream.h>
#endif

G4VVisCommandScene::G4VVisCommandScene () {}

void G4VVisCommandScene::UpdateCandidateLists () {

  G4SceneList& list = fpVisManager -> SetSceneDataObjectList ();
  // Seems it has to be non-const to satisfy iterator.
  G4SceneListIterator iterator (list);

  G4String fSceneNameList;
  while (++iterator) {
    fSceneNameList += iterator.key () + " ";
  }
  fSceneNameList = fSceneNameList.strip ();

  if (fpCommandSceneNotifyHandlers) {
    fpCommandSceneNotifyHandlers -> GetParameter (0) ->
      SetParameterCandidates (fSceneNameList);
  }

  if (fpCommandSceneRemove) {
    fpCommandSceneRemove -> GetParameter (0) ->
      SetParameterCandidates (fSceneNameList);
  }

  if (fpCommandSceneSelect) {
    fpCommandSceneSelect -> GetParameter (0) ->
      SetParameterCandidates (fSceneNameList);
  }

  if (fpCommandSceneHandlerAttach) {
    fpCommandSceneHandlerAttach -> GetParameter (0) ->
      SetParameterCandidates (fSceneNameList);
  }
}

void G4VVisCommandScene::UpdateVisManagerSceneDataAndViewParameters
(const G4String& sceneName) {
  G4Scene& currentScene = fpVisManager -> SetCurrentSceneData ();
  G4SceneList& sceneList = fpVisManager -> SetSceneDataObjectList ();
  if (sceneList.contains (sceneName)) {
    currentScene = sceneList [sceneName];
    G4ViewParameters& currentVP = fpVisManager -> SetCurrentViewParameters ();
    currentVP.SetCurrentTargetPoint (currentScene.GetStandardTargetPoint ());
    currentVP.SetZoomFactor (1.);
    currentVP.SetDolly (0.);
  }
  else {
    currentScene = G4Scene ("invalid scene");
  }
}

G4String G4VVisCommandScene::fSceneNameList;

////////////// /vis/scene/create ///////////////////////////////////////

G4VisCommandSceneCreate::G4VisCommandSceneCreate (): fId (0) {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/scene/create", this);
  fpCommand -> SetGuidance ("/vis/scene/create [<scene-name>]");
  fpCommand -> SetGuidance ("Creates an empty scene.");
  fpCommand -> SetGuidance ("Invents a name if not supplied.");
  fpCommand -> SetGuidance ("This scene becomes current.");
  fpCommand -> SetParameterName ("scene-name",
				 omitable = true,
				 currentAsDefault = true);
}

G4VisCommandSceneCreate::~G4VisCommandSceneCreate () {
  delete fpCommand;
}

G4String G4VisCommandSceneCreate::NextName () {
  char nextName [20];
  ostrstream ost (nextName, 20);
  ost << "scene-" << fId  << ends;
  return nextName;
}

G4String G4VisCommandSceneCreate::GetCurrentValue (G4UIcommand* command) {
  return NextName ();
}

void G4VisCommandSceneCreate::SetNewValue (G4UIcommand* command,
					   G4String newValue) {
  G4String& newName = newValue;
  G4String nextName = NextName ();

  if (newName == "") {
    newName = nextName;
  }
  if (newName == nextName) fId++;

  G4SceneList& list = fpVisManager -> SetSceneDataObjectList ();
  if (list.contains (newName)) {
    G4cout << "Scene \"" << newName << "\" already exists." << endl;
  }
  else {

    list [newName] = G4Scene (newName);
    // Adds empty scene data object to list.

    UpdateVisManagerSceneDataAndViewParameters (newName);
    G4cout << "New empty scene \"" << newName << "\" created." << endl;

    UpdateCandidateLists ();
  }
}

////////////// /vis/scene/list ///////////////////////////////////////

G4VisCommandSceneList::G4VisCommandSceneList () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/list", this);
  fpCommand -> SetGuidance ("/vis/scene/list [<scene-name>] [<verbosity>]");
  fpCommand -> SetGuidance ("Lists scene(s).");
  fpCommand -> SetGuidance ("<scene-name> default is \"all\"");
  fpCommand -> SetGuidance
    ("<verbosity> is 0 for short (default) or 1 for long listing.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("scene-name", 's',
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

G4VisCommandSceneList::~G4VisCommandSceneList () {
  delete fpCommand;
}

G4String G4VisCommandSceneList::GetCurrentValue (G4UIcommand* command) {
  return "";
}

void G4VisCommandSceneList::SetNewValue (G4UIcommand* command,
					 G4String newValue) {
  G4String name;
  G4int verbosity;
  istrstream is ((char*)newValue.data());
  is >> name >> verbosity;

  G4SceneList& list = fpVisManager -> SetSceneDataObjectList ();
  G4SceneListIterator iterator (list);
  G4bool found = false;
  while (++iterator) {
    if (name != "all") {
      if (name != iterator.key ()) continue;
    }
    found = true;
    G4cout << "Scene name: \"" << iterator.key () << "\"";
    if (verbosity > 0) {
      G4cout << ", " << iterator.value ();
    }
    G4cout << endl;
  }
  if (!found) {
    G4cout << "No scenes found";
    if (name != "all") {
      G4cout << " of name \"" << name << "\"";
    }
    G4cout << "." << endl;
  }
}

////////////// /vis/scene/notifyHandlers /////////////////////////

G4VisCommandSceneNotifyHandlers::G4VisCommandSceneNotifyHandlers () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/scene/notifyHandlers", this);
  fpCommand -> SetGuidance ("/vis/scene/notifyHandlers [<scene-name>]");
  fpCommand -> SetGuidance
    ("Notifies scene handlers of possible changes of scene.");
  fpCommand -> SetGuidance ("<scene-name> default is current scene.");
  fpCommand -> SetParameterName ("scene-name",
				 omitable = true,
				 currentAsDefault = true);
  fpCommandSceneNotifyHandlers = fpCommand;
}

G4VisCommandSceneNotifyHandlers::~G4VisCommandSceneNotifyHandlers () {
  delete fpCommand;
}

G4String G4VisCommandSceneNotifyHandlers::GetCurrentValue
(G4UIcommand* command) {
  return fpVisManager -> GetCurrentSceneData ().GetName ();
}

void G4VisCommandSceneNotifyHandlers::SetNewValue (G4UIcommand* command,
						   G4String newValue) {
  G4String& sceneName = newValue;
  const G4SceneList& sceneList =
    fpVisManager -> GetSceneDataObjectList ();
  G4SceneHandlerList& sceneHandlerList = fpVisManager -> SetAvailableScenes ();

  // For each scene handler, if it contains the scene, for each viewer
  // clear, (re)draw, and show.
  const G4int nSceneHandlers = sceneHandlerList.entries ();
  for (G4int iSH = 0; iSH < nSceneHandlers; iSH++) {
    G4VSceneHandler* aSceneHandler = sceneHandlerList [iSH];
    const G4Scene& theScene = aSceneHandler -> GetSceneData ();
    const G4String& theSceneName = theScene.GetName ();
    if (sceneName == theSceneName) {
      G4ViewerList& viewerList = aSceneHandler -> SetViewList ();
      const G4int nViewers = viewerList.entries ();
      for (G4int iV = 0; iV < nViewers; iV++) {
	G4VViewer* aViewer = viewerList [iV];
	aViewer -> ClearView ();
	aViewer -> DrawView ();
	aViewer -> ShowView ();
	G4cout << "Viewer \"" << aViewer -> GetName ()
	       << "\" of scene handler \"" << aSceneHandler -> GetName ()
	       << "\"\n  refreshed at request of scene \"" << sceneName
	       << "\"." << endl;
      }
    }
  }
}

////////////// /vis/scene/remove ///////////////////////////////////////

G4VisCommandSceneRemove::G4VisCommandSceneRemove () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/scene/remove", this);
  fpCommand -> SetGuidance ("/vis/scene/remove <scene-name>");
  fpCommand -> SetGuidance ("Removes scenes.");
  fpCommand -> SetGuidance
    ("/vis/scene/list to see possibilities.");
  fpCommand -> SetParameterName ("scene-name",
				 omitable = false,
				 currentAsDefault = true);
  fpCommandSceneRemove = fpCommand;
}

G4VisCommandSceneRemove::~G4VisCommandSceneRemove () {
  delete fpCommand;
}

G4String G4VisCommandSceneRemove::GetCurrentValue (G4UIcommand* command) {
  return fpVisManager -> GetCurrentSceneData ().GetName ();
}

void G4VisCommandSceneRemove::SetNewValue (G4UIcommand* command,
					   G4String newValue) {
  G4String& removeName = newValue;
  G4SceneList& list = fpVisManager -> SetSceneDataObjectList ();
  const G4String& currentSceneName =
    fpVisManager -> GetCurrentSceneData ().GetName ();
  if (list.contains (removeName)) {
    list.remove (removeName);
    G4cout << "Scene \"" << removeName << "\" removed." << endl;
    UpdateCandidateLists ();
    if (list.isEmpty ()) {
      UpdateVisManagerSceneDataAndViewParameters ();
      G4cout << "No scenes left.  Please create a new one." << endl;
    }
    else {
      if (currentSceneName == removeName) {
	G4SceneListIterator iterator (list);
	++iterator;
	UpdateVisManagerSceneDataAndViewParameters (iterator.key ());
	G4cout << "Current scene is now \""
	       << fpVisManager -> GetCurrentSceneData ().GetName ()
	       << "\"" << endl;
      }
      else {
	G4cout << "Current scene unchanged." << endl;
      }
    }
  }
  else {
    G4cout << "Scene \"" << removeName
	   << "\" not found - \"/vis/scene/list\" to see possibilities."
	   << endl;
  }
}

////////////// /vis/scene/select ///////////////////////////////////////

G4VisCommandSceneSelect::G4VisCommandSceneSelect () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/scene/select", this);
  fpCommand -> SetGuidance ("/vis/scene/select [<scene-name>]");
  fpCommand -> SetGuidance ("Selects current scene.");
  fpCommand -> SetGuidance
    ("Specify scene by name (\"/vis/scene/list\" to see possibilities).");
  fpCommand -> SetParameterName ("scene-name",
				 omitable = true,
				 currentAsDefault = true);
  fpCommandSceneSelect = fpCommand;
}

G4VisCommandSceneSelect::~G4VisCommandSceneSelect () {
  delete fpCommand;
}

G4String G4VisCommandSceneSelect::GetCurrentValue (G4UIcommand* command) {
  return fpVisManager -> GetCurrentSceneData ().GetName ();
}

void G4VisCommandSceneSelect::SetNewValue (G4UIcommand* command,
					   G4String newValue) {
  G4String& selectName = newValue;
  G4SceneList& list = fpVisManager -> SetSceneDataObjectList ();
  G4cout << "Scene \"" << selectName;
  if (list.contains (selectName)) {
    UpdateVisManagerSceneDataAndViewParameters (selectName);
    G4cout << "\" selected." << endl;
  }
  else {
    G4cout << "\" not found - \"/vis/scene/list\" to see possibilities."
	   << endl;
  }
}
