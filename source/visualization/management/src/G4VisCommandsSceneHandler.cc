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

// /vis/sceneHandler commands - John Allison  10th October 1998

#include "G4VisCommandsSceneHandler.hh"

#include "G4VisManager.hh"
#include "G4GraphicsSystemList.hh"
#include "G4VisCommandsScene.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"
#include <sstream>

#define G4warn G4cout

////////////// /vis/sceneHandler/attach ///////////////////////////////////////

G4VisCommandSceneHandlerAttach::G4VisCommandSceneHandlerAttach () {
  G4bool omitable, currentAsDefault;
  fpCommand = new G4UIcmdWithAString ("/vis/sceneHandler/attach", this);
  fpCommand -> SetGuidance ("Attaches scene to current scene handler.");
  fpCommand -> SetGuidance
    ("If scene-name is omitted, current scene is attached.  To see scenes and"
  "\nscene handlers, use \"/vis/scene/list\" and \"/vis/sceneHandler/list\"");
  fpCommand -> SetParameterName ("scene-name",
				 omitable = true,
				 currentAsDefault = true);
}

G4VisCommandSceneHandlerAttach::~G4VisCommandSceneHandlerAttach () {
  delete fpCommand;
}

G4String G4VisCommandSceneHandlerAttach::GetCurrentValue (G4UIcommand*) {
  G4Scene* pScene = fpVisManager -> GetCurrentScene ();
  return pScene ? pScene -> GetName () : G4String("");
}

void G4VisCommandSceneHandlerAttach::SetNewValue (G4UIcommand*,
						  G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String& sceneName = newValue;

  if (sceneName.length () == 0) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout <<
	"WARNING: No scene specified.  Maybe there are no scenes available"
	"\n  yet.  Please create one." << G4endl;
    }
    return;
  }

  G4VSceneHandler* pSceneHandler = fpVisManager -> GetCurrentSceneHandler ();
  if (!pSceneHandler) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
      "ERROR: Current scene handler not defined.  Please select or create one."
	     << G4endl;
    }
    return;
  }

  G4SceneList& sceneList = fpVisManager -> SetSceneList ();

  if (sceneList.empty ()) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
      "ERROR: No valid scenes available yet.  Please create one."
	     << G4endl;
    }
    return;
  }

  std::size_t iScene, nScenes = sceneList.size ();
  for (iScene = 0; iScene < nScenes; ++iScene) {
    if (sceneList [iScene] -> GetName () == sceneName) break;
  }
  if (iScene < nScenes) {
    G4Scene* pScene = sceneList [iScene];
    pSceneHandler -> SetScene (pScene);
    // Make sure scene is current...
    fpVisManager -> SetCurrentScene (pScene);
    // Refresh viewer, if any (only if auto-refresh)...
    G4VViewer* pViewer = pSceneHandler -> GetCurrentViewer();
    if (pViewer && pViewer -> GetViewParameters().IsAutoRefresh()) {
      pViewer -> SetView ();
      pViewer -> ClearView ();
      pViewer -> DrawView ();
    }
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Scene \"" << sceneName
	     << "\" attached to scene handler \""
	     << pSceneHandler -> GetName () <<
	".\n  (You may have to refresh with \"/vis/viewer/flush\" if view"
	" is not \"auto-refresh\".)"
	     << G4endl;
    }
  }
  else {
    if (verbosity >= G4VisManager::errors) {
      G4warn << "ERROR: Scene \"" << sceneName
	     << "\" not found.  Use \"/vis/scene/list\" to see possibilities."
	     << G4endl;
    }
  }
}

////////////// /vis/sceneHandler/create ///////////////////////////////////////

G4VisCommandSceneHandlerCreate::G4VisCommandSceneHandlerCreate (): fId (0) {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/sceneHandler/create", this);
  fpCommand -> SetGuidance
    ("Creates an scene handler for a specific graphics system.");
  fpCommand -> SetGuidance
    ("Attaches current scene, if any.  (You can change attached scenes with"
     "\n\"/vis/sceneHandler/attach\".)  Invents a scene handler name if not"
     "\nsupplied.  This scene handler becomes current.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("graphics-system-name",
				 's', omitable = false);
  const G4GraphicsSystemList& gslist =
  fpVisManager -> GetAvailableGraphicsSystems ();
  G4String candidates;
  for (const auto gs: gslist) {
    const G4String& name = gs -> GetName ();
    candidates += name + ' ';
    for (const auto& nickname: gs -> GetNicknames ()) {
      if (G4StrUtil::contains(nickname, "FALLBACK")) continue;
      if (nickname != name) candidates += nickname + ' ';
    }
  }
  G4StrUtil::strip(candidates);
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
  std::ostringstream oss;
  oss << "scene-handler-" << fId;
  return oss.str();
}

G4String G4VisCommandSceneHandlerCreate::GetCurrentValue(G4UIcommand*) {

  G4String graphicsSystemName;
  const G4VGraphicsSystem* graphicsSystem =
    fpVisManager -> GetCurrentGraphicsSystem ();
  if (graphicsSystem) {
    graphicsSystemName = graphicsSystem -> GetName ();
  }
  else {
    const G4GraphicsSystemList& gslist =
      fpVisManager -> GetAvailableGraphicsSystems ();
    if (gslist.size ()) {
      graphicsSystemName = gslist [0] -> GetName ();
    }
    else {
      graphicsSystemName = "none";
    }
  }

  return graphicsSystemName + " " + NextName ();
}

void G4VisCommandSceneHandlerCreate::SetNewValue (G4UIcommand* command,
						  G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String graphicsSystem, newName;
  std::istringstream is (newValue);
  is >> graphicsSystem >> newName;

  const G4GraphicsSystemList& gsl =
    fpVisManager -> GetAvailableGraphicsSystems ();
  std::size_t nSystems = gsl.size ();
  if (nSystems <= 0) {
    G4ExceptionDescription ed;
    ed <<
    "ERROR: G4VisCommandSceneHandlerCreate::SetNewValue:"
    " no graphics systems available."
    "\n  Did you instantiate any in"
    " YourVisManager::RegisterGraphicsSystems()?";
    command->CommandFailed(ed);
    return;
  }
  std::size_t iGS;  // Selector index.
  G4bool found = false;
  for (iGS = 0; iGS < nSystems; ++iGS) {
    const auto& gs = gsl[iGS];
    if (G4StrUtil::icompare(graphicsSystem, gs->GetName()) == 0) {
      found = true;
      break;  // Match found
    } else {
      const auto& nicknames = gs->GetNicknames();
      for (std::size_t i = 0; i < nicknames.size(); ++i) {
        const auto& nickname = nicknames[i];
        if (G4StrUtil::icompare(graphicsSystem, nickname) == 0) {
          found = true;
          break;  // Match found
        }
      }
      if (found) {
        break;  // Match found
      }
    }
  }
  if (!found) {
    // Shouldn't happen, since graphicsSystem should be a candidate
    G4ExceptionDescription ed;
    ed <<
    "ERROR: G4VisCommandSceneHandlerCreate::SetNewValue:"
    "\n  Invalid graphics system \""
    << graphicsSystem
    << "\" requested."
    << "\n  Candidates are:";
    fpVisManager->PrintAvailableGraphicsSystems(verbosity,ed);
    command->CommandFailed(ed);
    return;
  }

  // Check UI session compatibility.
  G4bool fallback = false;
  G4int loopCounter = 0;
  while (!gsl[iGS]->IsUISessionCompatible()) {
    std::size_t iGSBeingTested = iGS;
    // Not compatible, search for a fallback
    fallback = false;
    G4String fallbackNickname = gsl[iGS]->GetNickname() + "_FALLBACK";
    for (iGS = 0; iGS < nSystems; iGS++) {
      const auto& nicknames = gsl[iGS]->GetNicknames();
      for (std::size_t i = 0; i < nicknames.size(); ++i) {
        const auto& nickname = nicknames[i];
        if (G4StrUtil::icompare(fallbackNickname, nickname) == 0) {
          fallback = true;
          break;  // Match found
        }
      }
      if (fallback) {
        break;  // Match found
      }
    }
    if (iGS >= nSystems || loopCounter >=3) {
      G4ExceptionDescription ed;
      ed << "\"" << gsl[iGSBeingTested]->GetNickname()
      << "\" is not compatible with your chosen session,"
      " and no fallback system found.";
      command->CommandFailed(ed);
      return;
    }
    //  A fallback system found...but go back and check this too.
    ++loopCounter;
  }

  // A graphics system has been found
  G4VGraphicsSystem* pSystem = gsl [iGS];

  if (fallback && verbosity >= G4VisManager::warnings) {
    G4warn << "WARNING: G4VisCommandSceneHandlerCreate::SetNewValue:"
    "\n  Using fallback graphics system: "
    << pSystem -> GetName ()
    << " ("
    << pSystem -> GetNickname ()
    << ')'
    << G4endl;
  }

  // Now deal with name of scene handler.
  G4String nextName = NextName ();
  if (newName == "") {
    newName = nextName;
  }
  if (newName == nextName) fId++;

  const G4SceneHandlerList& list = fpVisManager -> GetAvailableSceneHandlers ();
  std::size_t iScene;
  for (iScene = 0; iScene < list.size (); ++iScene) {
    G4VSceneHandler* sceneHandler = list [iScene];
    if (sceneHandler -> GetName () == newName) {
      G4ExceptionDescription ed;
      ed <<
      "ERROR: Scene handler \"" << newName
      << "\" already exists.";
      command->CommandFailed(ed);
      return;
    }
  }

  // If there is an existing viewer, store its view parameters
  if (fpVisManager->GetCurrentViewer()) {
    fThereWasAViewer = true;
    fExistingVP = fpVisManager->GetCurrentViewer()->GetViewParameters();
  }

  // Set current graphics system in preparation for
  // creating scene handler.
  fpVisManager -> SetCurrentGraphicsSystem (pSystem);
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Graphics system set to "
    << pSystem -> GetName ()
    << " ("
    << pSystem -> GetNickname ()
    << ')'
    << G4endl;
  }

  //Create scene handler.
  fpVisManager -> CreateSceneHandler (newName);
  if (fpVisManager -> GetCurrentSceneHandler () -> GetName () != newName) {
    G4ExceptionDescription ed;
    ed <<
    "ERROR: G4VisCommandSceneHandlerCreate::SetNewValue:"
    " Curious name mismatch."
    "\n Current name \""
    << fpVisManager -> GetCurrentSceneHandler () -> GetName ()
    << "\" is not the new name \""
    << newName
    << "\".\n  Please report to vis coordinator.";
    command->CommandFailed(ed);
    return;
  }

  if (verbosity >= G4VisManager::confirmations)
    G4cout << "New scene handler \"" << newName << "\" created." << G4endl;

  if (fpVisManager -> GetCurrentScene ()) {
    auto errorCode = G4UImanager::GetUIpointer () -> ApplyCommand ("/vis/sceneHandler/attach");
    if (errorCode) {
      G4ExceptionDescription ed;
      ed << "sub-command \"/vis/sceneHandler/attach\" failed.";
      command->CommandFailed(errorCode,ed);
      return;
    }
  }
}

////////////// /vis/sceneHandler/list ///////////////////////////////////////

G4VisCommandSceneHandlerList::G4VisCommandSceneHandlerList () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/sceneHandler/list", this);
  fpCommand -> SetGuidance ("Lists scene handler(s).");
  fpCommand -> SetGuidance
    ("\"help /vis/verbose\" for definition of verbosity.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter("scene-handler-name", 's', omitable = true);
  parameter -> SetDefaultValue ("all");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("verbosity", 's', omitable = true);
  parameter -> SetDefaultValue ("warnings");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneHandlerList::~G4VisCommandSceneHandlerList () {
  delete fpCommand;
}

G4String G4VisCommandSceneHandlerList::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneHandlerList::SetNewValue (G4UIcommand*,
						G4String newValue) {
  G4String name, verbosityString;
  std::istringstream is (newValue);
  is >> name >> verbosityString;
  G4VisManager::Verbosity verbosity =
    fpVisManager->GetVerbosityValue(verbosityString);
  const G4VSceneHandler* currentSceneHandler =
    fpVisManager -> GetCurrentSceneHandler ();
  G4String currentName;
  if (currentSceneHandler) currentName = currentSceneHandler->GetName();

  const G4SceneHandlerList& list = fpVisManager -> GetAvailableSceneHandlers ();
  G4bool found = false;
  for (std::size_t iSH = 0; iSH < list.size (); ++iSH) {
    const G4String& iName = list [iSH] -> GetName ();
    if (name != "all") {
      if (name != iName) continue;
    }
    found = true;
    if (iName == currentName) {
      G4cout << "  (current)";
    }
    else {
      G4cout << "           ";
    }
    G4cout << " scene handler \"" << list [iSH] -> GetName () << "\""
	   << " (" << list [iSH] -> GetGraphicsSystem () -> GetName () << ")";
    if (verbosity >= G4VisManager::parameters) {
      G4cout << "\n  " << *(list [iSH]);
    }
    G4cout << G4endl;
  }
  if (!found) {
    G4cout << "No scene handlers found";
    if (name != "all") {
      G4cout << " of name \"" << name << "\"";
    }
    G4cout << "." << G4endl;
  }
}

////////////// /vis/sceneHandler/select ///////////////////////////////////////

G4VisCommandSceneHandlerSelect::G4VisCommandSceneHandlerSelect () {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString ("/vis/sceneHandler/select", this);
  fpCommand -> SetGuidance ("Selects a scene handler.");
  fpCommand -> SetGuidance 
    ("Makes the scene handler current.  \"/vis/sceneHandler/list\" to see"
     "\n possible scene handler names.");
  fpCommand -> SetParameterName ("scene-handler-name",
				 omitable = false);
}

G4VisCommandSceneHandlerSelect::~G4VisCommandSceneHandlerSelect () {
  delete fpCommand;
}

G4String G4VisCommandSceneHandlerSelect::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneHandlerSelect::SetNewValue (G4UIcommand*,
						  G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String& selectName = newValue;
  const G4SceneHandlerList& list = fpVisManager -> GetAvailableSceneHandlers ();

  std::size_t iSH;
  for (iSH = 0; iSH < list.size (); iSH++) {
    if (list [iSH] -> GetName () == selectName) break;
  }
  if (iSH < list.size ()) {
    if (fpVisManager -> GetCurrentSceneHandler () -> GetName ()
	== selectName) {
      if (verbosity >= G4VisManager::confirmations) {
	G4cout << "Scene handler \"" << selectName << "\""
	       << " already selected." << G4endl;
      }
    }
    else {
      if (verbosity >= G4VisManager::confirmations) {
	G4cout << "Scene handler \"" << selectName << "\""
	       << " being selected." << G4endl;
      }
      fpVisManager -> SetCurrentSceneHandler (list [iSH]);
    }
  }
  else {
    if (verbosity >= G4VisManager::errors) {
      G4warn << "ERROR: Scene handler \"" << selectName << "\""
	     << " not found - \"/vis/sceneHandler/list\" to see possibilities."
	     << G4endl;
    }
  }
}
