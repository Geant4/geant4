// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsSceneHandler.cc,v 1.6 1999-05-10 14:04:17 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/sceneHandler commands - John Allison  10th October 1998

#include "G4VisCommandsSceneHandler.hh"

#include "G4VisManager.hh"
#include "G4GraphicsSystemList.hh"
#include "G4VisCommandsScene.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"
#ifdef WIN32
#include <strstrea.h>
#else
#include <strstream.h>
#endif

G4VVisCommandSceneHandler::G4VVisCommandSceneHandler () {}

G4VVisCommandSceneHandler::~G4VVisCommandSceneHandler () {}

void G4VVisCommandSceneHandler::UpdateCandidateLists () {

  const G4SceneHandlerList& list =
    fpVisManager -> GetAvailableSceneHandlers ();

  fSceneHandlerNameList = G4String ();
  for (int iScene = 0; iScene < list.entries (); iScene++) {
    G4VSceneHandler* sceneHandler = list [iScene];
    fSceneHandlerNameList += sceneHandler -> GetName () + " ";
  }
  fSceneHandlerNameList = fSceneHandlerNameList.strip ();

  if (fpCommandSceneHandlerRemove) {
    fpCommandSceneHandlerRemove -> GetParameter (0) ->
      SetParameterCandidates (fSceneHandlerNameList);
  }

  if (fpCommandSceneHandlerSelect) {
    fpCommandSceneHandlerSelect -> GetParameter (0) ->
      SetParameterCandidates (fSceneHandlerNameList);
  }

  if (fpCommandViewerCreate) {
    fpCommandViewerCreate -> GetParameter (0) ->
      SetParameterCandidates (fSceneHandlerNameList);
  }
}

G4String G4VVisCommandSceneHandler::fSceneHandlerNameList;

////////////// /vis/sceneHandler/attach ///////////////////////////////////////

G4VisCommandSceneHandlerAttach::G4VisCommandSceneHandlerAttach () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/sceneHandler/attach", this);
  fpCommand -> SetGuidance
    ("/vis/sceneHandler/attach [<scene-name>]");
  fpCommand -> SetGuidance ("Attaches scene to current scene handler.");
  fpCommand -> SetGuidance
    ("If scene-name is omitted, current scene is attached.");
  fpCommand -> SetGuidance
    ("To see scenes and scene handlers, use \"/vis/scene/list\""
     "\n  and \"/vis/sceneHandler/list\"");
  fpCommand -> SetParameterName ("scene-name",
				 omitable = true,
				 currentAsDefault = true);
  fpCommandSceneHandlerAttach = fpCommand;
}

G4VisCommandSceneHandlerAttach::~G4VisCommandSceneHandlerAttach () {
  delete fpCommand;
}

G4String G4VisCommandSceneHandlerAttach::GetCurrentValue (G4UIcommand* command) {
  return fpVisManager -> GetCurrentScene () -> GetName ();
}

void G4VisCommandSceneHandlerAttach::SetNewValue (G4UIcommand* command,
						  G4String newValue) {
  G4String& sceneName = newValue;

  if (sceneName.length () == 0) {
    G4cout <<
      "Null string specified.  Maybe there are no scenes available yet."
      "\n  Please create one." << endl;
    return;
  }

  G4VSceneHandler* pSceneHandler = fpVisManager -> GetCurrentSceneHandler ();
  if (!pSceneHandler) {
    G4cout <<
      "Current scene handler not defined.  Please select or create one."
	   << endl;
    return;
  }

  G4SceneList& sceneList =
    fpVisManager -> SetSceneList ();

  if (sceneList.isEmpty ()) {
    G4cout <<
      "No valid scenes available yet.  Please create one."
	   << endl;
    return;
  }

  G4int iScene, nScenes = sceneList.entries ();
  for (iScene = 0; iScene < nScenes; iScene++) {
    if (sceneList [iScene] -> GetName () == sceneName) break;
  }
  if (iScene < nScenes) {
    G4Scene* pScene = sceneList [iScene];
    pSceneHandler -> SetScene (pScene);
    fpVisManager -> SetCurrentScene (pScene);
    G4cout << "Scene \"" << sceneName
	   << "\" attached to scene handler \"" << pSceneHandler -> GetName ()
	   << "." << endl;
  }
  else {
    G4cout << "Scene \"" << sceneName
	   << "\" not found.  Use \"/vis/scene/list\" to see possibilities."
	   << endl;
  }
}

////////////// /vis/sceneHandler/create ///////////////////////////////////////

G4VisCommandSceneHandlerCreate::G4VisCommandSceneHandlerCreate (): fId (0) {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/sceneHandler/create", this);
  fpCommand -> SetGuidance
    ("/vis/sceneHandler/create");
  fpCommand -> SetGuidance
    ("     <graphics-system> [<scene-handler-name>] [<scene-name>]");
  fpCommand -> SetGuidance
    ("Creates an scene handler for a specific graphics system.");
  fpCommand -> SetGuidance
    ("Attaches specified scene.");
  fpCommand -> SetGuidance
    ("Default graphics system is current graphics system.");
  fpCommand -> SetGuidance
    ("Invents a name if not supplied.");
  fpCommand -> SetGuidance
    ("Default scene is current scene.  You can change attached scenes with");
  fpCommand -> SetGuidance
    ("  /vis/sceneHandler/attach [<scene-name>].");
  fpCommand -> SetGuidance
    ("This graphics system, scene handler and scene become current.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("graphics-system", 's', omitable = false);
  //parameter -> SetCurrentAsDefault (true);
  const G4GraphicsSystemList& gslist =
    fpVisManager -> GetAvailableGraphicsSystems ();
  G4String candidates;
  for (int igslist = 0; igslist < gslist.entries (); igslist++) {
    const G4String& name = gslist (igslist) -> GetName ();
    const G4String& nickname = gslist (igslist) -> GetNickname ();
    if (nickname.isNull ()) {
      candidates += name;
    }
    else {
      candidates += nickname;
    }
    candidates += " ";
  }
  candidates = candidates.strip ();
  parameter -> SetParameterCandidates(candidates);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter
    ("scene-handler-name", 's', omitable = true);
  parameter -> SetCurrentAsDefault (true);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneHandlerCreate::~G4VisCommandSceneHandlerCreate () {
  delete fpCommand;
}

G4String G4VisCommandSceneHandlerCreate::NextName () {
  char nextName [20];
  ostrstream ost (nextName, 20);
  ost << "scene-handler-" << fId << ends;
  return nextName;
}

G4String G4VisCommandSceneHandlerCreate::GetCurrentValue
(G4UIcommand* command) {

  G4String graphicsSystemName;
  const G4VGraphicsSystem* graphicsSystem =
    fpVisManager -> GetCurrentGraphicsSystem ();
  if (graphicsSystem) {
    graphicsSystemName = graphicsSystem -> GetName ();
  }
  else {
    const G4GraphicsSystemList& gslist =
      fpVisManager -> GetAvailableGraphicsSystems ();
    if (gslist.entries ()) {
      graphicsSystemName = gslist [0] -> GetName ();
    }
    else {
      graphicsSystemName = "none";
    }
  }

  G4String sceneName;
  const G4VSceneHandler* scene = fpVisManager -> GetCurrentSceneHandler ();
  if (scene) {
    sceneName = scene -> GetName ();
  }
  else {
    sceneName = "none";
  }

  return graphicsSystemName + " " + NextName () + " " + sceneName;
}

void G4VisCommandSceneHandlerCreate::SetNewValue (G4UIcommand* command,
					   G4String newValue) {
  G4String graphicsSystem, newName, sceneName;
  istrstream is ((char*)newValue.data());
  is >> graphicsSystem >> newName >> sceneName;

  const G4GraphicsSystemList& gsl =
    fpVisManager -> GetAvailableGraphicsSystems ();
  int nSystems = gsl.entries ();
  if (nSystems <= 0) {
    G4cout << "G4VisCommandSceneHandlerCreate::SetNewValue:"
      " no graphics systems available."
      "\n  Did you instantiate any in"
      " YourVisManager::RegisterGraphicsSystems()?"
	   << endl;
    return;
  }
  int iGS;  // Selector index.
  for (iGS = 0; iGS < nSystems; iGS++) {
    if (graphicsSystem.compareTo (gsl [iGS] -> GetName (),
				  RWCString::ignoreCase) == 0 ||
	graphicsSystem.compareTo (gsl [iGS] -> GetNickname (),
				  RWCString::ignoreCase) == 0) {
      break;  // Match found.
    }
  }
  if (iGS < 0 || iGS >= nSystems) {
    // Invalid command line argument or non.
    // This shouldn't happen!!!!!!
    G4cerr << "G4VisCommandSceneHandlerCreate::SetNewValue:"
      " invalid graphics system specified."
	   << endl;
    return;
  }
  // Valid index.  Set current graphics system in preparation for
  // creating scene handler.
  G4VGraphicsSystem* pSystem = gsl [iGS];
  fpVisManager -> SetCurrentGraphicsSystem (pSystem);
  if (fpVisManager -> GetVerboseLevel () > 0) {
    G4cout << "Graphics system set to " << pSystem -> GetName () << endl;
  }

  // Now deal with name of scene handler.
  G4String nextName = NextName ();
  if (newName == "") {
    newName = nextName;
  }
  if (newName == nextName) fId++;

  const G4SceneHandlerList& list = fpVisManager -> GetAvailableSceneHandlers ();
  int iScene;
  for (iScene = 0; iScene < list.entries (); iScene++) {
    G4VSceneHandler* sceneHandler = list [iScene];
    if (sceneHandler -> GetName () == newName) {
      G4cout << "Scene handler \"" << newName << "\" already exists." << endl;
      return;
    }
  }

  //Create scene handler.
  fpVisManager -> CreateSceneHandler (newName);
  G4cout << "New scene handler \"" << newName << "\" created." << endl;

  // Attach scene.
  G4String commandString = "/vis/sceneHandler/attach " + sceneName;
  G4UImanager::GetUIpointer () -> ApplyCommand (commandString);

  UpdateCandidateLists ();
}

////////////// /vis/sceneHandler/list ///////////////////////////////////////

G4VisCommandSceneHandlerList::G4VisCommandSceneHandlerList () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/sceneHandler/list", this);
  fpCommand -> SetGuidance
    ("/vis/sceneHandler/list [<scene-handler-name>] [<verbosity>]");
  fpCommand -> SetGuidance ("Lists scene handler(s).");
  fpCommand -> SetGuidance ("<scene-handler-name> default is \"all\"");
  fpCommand -> SetGuidance
    ("<verbosity> is 0 for short (default) or 1 for long listing.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter("scene-handler-name", 's',
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

G4VisCommandSceneHandlerList::~G4VisCommandSceneHandlerList () {
  delete fpCommand;
}

G4String G4VisCommandSceneHandlerList::GetCurrentValue (G4UIcommand* command) {
  return "";
}

void G4VisCommandSceneHandlerList::SetNewValue (G4UIcommand* command,
					 G4String newValue) {
  G4String name;
  G4int verbosity;
  istrstream is ((char*)newValue.data());
  is >> name >> verbosity;

  const G4SceneHandlerList& list = fpVisManager -> GetAvailableSceneHandlers ();
  G4bool found = false;
  for (int iSH = 0; iSH < list.entries (); iSH++) {
    if (name != "all") {
      if (name != list [iSH] -> GetName ()) continue;
    }
    found = true;
    G4cout << "Scene handler \"" << list [iSH] -> GetName () << "\""
	   << " (" << list [iSH] -> GetGraphicsSystem () -> GetName () << ")";
    if (verbosity > 0) {
      G4cout << "\n  " << *(list [iSH]);
    }
    G4cout << endl;
  }
  if (!found) {
    G4cout << "No scene handlers found";
    if (name != "all") {
      G4cout << " of name \"" << name << "\"";
    }
    G4cout << "." << endl;
  }
}

////////////// /vis/sceneHandler/remove ///////////////////////////////////////

G4VisCommandSceneHandlerRemove::G4VisCommandSceneHandlerRemove () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/sceneHandler/remove", this);
  fpCommand -> SetGuidance ("/vis/sceneHandler/remove <scene-handler-name>");
  fpCommand -> SetGuidance ("Removes scene handlers.");
  fpCommand -> SetGuidance
    ("Specify scene handler by name (\"/vis/sceneHandler/list\""
     "\n  to see possibilities).");
  fpCommand -> SetParameterName ("scene-handler-name",
				 omitable = false,
				 currentAsDefault = true);
  fpCommandSceneHandlerRemove = fpCommand;
}

G4VisCommandSceneHandlerRemove::~G4VisCommandSceneHandlerRemove () {
  delete fpCommand;
}

G4String G4VisCommandSceneHandlerRemove::GetCurrentValue (G4UIcommand* command) {
  G4VSceneHandler* sceneHandler = fpVisManager -> GetCurrentSceneHandler ();
  if (sceneHandler) {
    return sceneHandler -> GetName ();
  }
  else {
    return "none";
  }
}

void G4VisCommandSceneHandlerRemove::SetNewValue (G4UIcommand* command,
					   G4String newValue) {
  G4String& removeName = newValue;
  G4VSceneHandler* currentSceneHandler =
    fpVisManager -> GetCurrentSceneHandler ();
  G4String currentName;
  if (currentSceneHandler) {
    currentName = currentSceneHandler -> GetName ();
  }

  G4SceneHandlerList& list = fpVisManager -> SetAvailableSceneHandlers ();
  G4int iSH;
  for (iSH = 0; iSH < list.entries (); iSH++) {
    if (list [iSH] -> GetName () == removeName) break;
  }

  if (iSH >= list.entries ()) {
    G4cout << "Scene handler \"" << removeName
	   << "\" not found - \"/vis/sceneHandler/list\" to see possibilities."
	   << endl;
    return;
  }

  G4cout << "Scene handler \"" << removeName << "\" removed." << endl;
  if (removeName == currentName) {
    fpVisManager -> DeleteCurrentSceneHandler ();
  }
  else {
    list.remove (list [iSH]);
    G4cout << "Current scene handler unchanged." << endl;
  }

  UpdateCandidateLists ();
}

////////////// /vis/sceneHandler/select ///////////////////////////////////////

G4VisCommandSceneHandlerSelect::G4VisCommandSceneHandlerSelect () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/sceneHandler/select", this);
  fpCommand -> SetGuidance ("/vis/sceneHandler/select [<scene-handler-name>]");
  fpCommand -> SetGuidance ("Selects current scene handler.");
  fpCommand -> SetGuidance
    ("Specify scene handler by name (\"/vis/sceneHandler/list\""
     "\n  to see possibilities).");
  fpCommand -> SetParameterName ("scene-handler-name",
				 omitable = true,
				 currentAsDefault = true);
  fpCommandSceneHandlerSelect = fpCommand;
}

G4VisCommandSceneHandlerSelect::~G4VisCommandSceneHandlerSelect () {
  delete fpCommand;
}

G4String G4VisCommandSceneHandlerSelect::GetCurrentValue 
(G4UIcommand* command) {
  G4VSceneHandler* sceneHandler = fpVisManager -> GetCurrentSceneHandler ();
  if (sceneHandler) {
    return sceneHandler -> GetName ();
  }
  else {
    return "none";
  }
}

void G4VisCommandSceneHandlerSelect::SetNewValue (G4UIcommand* command,
					   G4String newValue) {
  G4String& selectName = newValue;
  const G4SceneHandlerList& list = fpVisManager -> GetAvailableSceneHandlers ();
  G4cout << "Scene handler \"" << selectName << "\"";
  G4int iSH;
  for (iSH = 0; iSH < list.entries (); iSH++) {
    if (list [iSH] -> GetName () == selectName) break;
  }
  if (iSH < list.entries ()) {
    if (fpVisManager -> GetCurrentSceneHandler () -> GetName () == selectName) {
      G4cout << " already selected." << endl;
    }
    else {
      G4cout << " being selected." << endl;
      fpVisManager -> SetCurrentSceneHandler (list [iSH]);
    }
  }
  else {
    G4cout << " not found - \"/vis/sceneHandler/list\""
      "\n  to see possibilities."
	   << endl;
  }
}
