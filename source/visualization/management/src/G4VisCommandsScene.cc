// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsScene.cc,v 1.5 1999-03-29 16:39:11 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/scene commands - John Allison  9th August 1998

#include "G4VisCommandsScene.hh"

#include "G4VisManager.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4ApplicationState.hh"
#include "G4UImanager.hh"
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

  G4SceneList& sceneList = fpVisManager -> SetSceneList ();
  G4int iScene, nScenes = sceneList.entries ();
  G4String nameList;
  for (iScene = 0; iScene < nScenes; iScene++) {
    nameList += sceneList [iScene] -> GetName () + " ";
  }
  nameList = nameList.strip ();

  if (fpCommandSceneNotifyHandlers) {
    fpCommandSceneNotifyHandlers -> GetParameter (0) ->
      SetParameterCandidates (nameList);
  }

  if (fpCommandSceneRemove) {
    fpCommandSceneRemove -> GetParameter (0) ->
      SetParameterCandidates (nameList);
  }

  if (fpCommandSceneSelect) {
    fpCommandSceneSelect -> GetParameter (0) ->
      SetParameterCandidates (nameList);
  }

  if (fpCommandSceneHandlerAttach) {
    fpCommandSceneHandlerAttach -> GetParameter (0) ->
      SetParameterCandidates (nameList);
  }
}

void G4VVisCommandScene::UpdateVisManagerSceneAndViewParameters
(const G4String& sceneName) {
  G4SceneList& sceneList = fpVisManager -> SetSceneList ();
  G4int iScene, nScenes = sceneList.entries ();
  for (iScene = 0; iScene < nScenes; iScene++) {
    if (sceneList [iScene] -> GetName () == sceneName) break;
  }
  if (iScene < nScenes) {
    fpVisManager -> SetCurrentScene (sceneList [iScene]);
    G4UImanager::GetUIpointer () -> ApplyCommand ("/vis/camera/reset");
  }
  else {
    fpVisManager -> SetCurrentScene (0);
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

  G4SceneList& sceneList = fpVisManager -> SetSceneList ();
  G4int iScene, nScenes = sceneList.entries ();
  for (iScene = 0; iScene < nScenes; iScene++) {
    if (sceneList [iScene] -> GetName () == newName) break;
  }
  if (iScene < nScenes) {
    G4cout << "Scene \"" << newName << "\" already exists." << endl;
  }
  else {

    sceneList.append (new G4Scene (newName));
    // Adds empty scene data object to list.

    UpdateVisManagerSceneAndViewParameters (newName);
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

  G4SceneList& sceneList = fpVisManager -> SetSceneList ();
  G4int iScene, nScenes = sceneList.entries ();
  G4bool found = false;
  for (iScene = 0; iScene < nScenes; iScene++) {
    const G4String& iName = sceneList [iScene] -> GetName ();
    if (name != "all") {
      if (name != iName) continue;
    }
    found = true;
    G4cout << "Scene name: \"" << iName << "\"";
    if (verbosity > 0) {
      G4cout << ", " << *sceneList [iScene];
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
  fpCommand -> SetGuidance ("<scene-name> default is current scene name.");
  fpCommand -> SetGuidance
    ("This command does not change current scene, scene handler or viewer.");
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
  return fpVisManager -> GetCurrentScene () -> GetName ();
}

void G4VisCommandSceneNotifyHandlers::SetNewValue (G4UIcommand* command,
						   G4String newValue) {
  G4String& sceneName = newValue;
  const G4SceneList& sceneList = fpVisManager -> GetSceneList ();
  G4SceneHandlerList& sceneHandlerList =
    fpVisManager -> SetAvailableSceneHandlers ();

  // Check scene name.
  const G4int nScenes = sceneList.entries ();
  G4int iScene;
  for (iScene = 0; iScene < nScenes; iScene++) {
    G4Scene* scene = sceneList [iScene];
    if (sceneName == scene -> GetName ()) {
      scene -> AddWorldIfEmpty ();
      break;
    }
  }
  if (iScene >= nScenes ) {
    G4cout << "Scene \"" << sceneName << "\" not found."
      "\n  /vis/scene/list to see scenes."
	   << endl;
    return;
  }

  // For each scene handler, if it contains the scene, for each viewer
  // set (make current), clear, (re)draw, and show.
  const G4int nSceneHandlers = sceneHandlerList.entries ();
  for (G4int iSH = 0; iSH < nSceneHandlers; iSH++) {
    G4VSceneHandler* aSceneHandler = sceneHandlerList [iSH];
    G4Scene* theScene = aSceneHandler -> GetScene ();
    const G4String& theSceneName = theScene -> GetName ();
    if (sceneName == theSceneName) {
      G4ViewerList& viewerList = aSceneHandler -> SetViewerList ();
      const G4int nViewers = viewerList.entries ();
      for (G4int iV = 0; iV < nViewers; iV++) {
	G4VViewer* aViewer = viewerList [iV];
	aViewer -> SetView ();
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
  return fpVisManager -> GetCurrentScene () -> GetName ();
}

void G4VisCommandSceneRemove::SetNewValue (G4UIcommand* command,
					   G4String newValue) {
  G4String& removeName = newValue;
  const G4String& currentSceneName =
    fpVisManager -> GetCurrentScene () -> GetName ();
  G4SceneList& sceneList = fpVisManager -> SetSceneList ();
  G4int iScene, nScenes = sceneList.entries ();
  for (iScene = 0; iScene < nScenes; iScene++) {
    if (sceneList [iScene] -> GetName () == removeName) break;
  }
  if (iScene < nScenes) {
    delete sceneList [iScene];
    sceneList.removeAt (iScene);
    G4cout << "Scene \"" << removeName << "\" removed." << endl;
    UpdateCandidateLists ();
    if (sceneList.isEmpty ()) {
      UpdateVisManagerSceneAndViewParameters ();
      G4cout << "No scenes left.  Please create a new one." << endl;
    }
    else {
      if (currentSceneName == removeName) {
	UpdateVisManagerSceneAndViewParameters (sceneList [0] -> GetName ());
	G4cout << "Current scene is now \""
	       << fpVisManager -> GetCurrentScene () -> GetName ()
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
  return fpVisManager -> GetCurrentScene () -> GetName ();
}

void G4VisCommandSceneSelect::SetNewValue (G4UIcommand* command,
					   G4String newValue) {
  G4String& selectName = newValue;
  G4SceneList& sceneList = fpVisManager -> SetSceneList ();
  G4int iScene, nScenes = sceneList.entries ();
  for (iScene = 0; iScene < nScenes; iScene++) {
    if (sceneList [iScene] -> GetName () == selectName) break;
  }
  G4cout << "Scene \"" << selectName;
  if (iScene < nScenes) {
    UpdateVisManagerSceneAndViewParameters (selectName);
    G4cout << "\" selected." << endl;
  }
  else {
    G4cout << "\" not found - \"/vis/scene/list\" to see possibilities."
	   << endl;
  }
}
