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

// /vis/scene commands - John Allison  9th August 1998

#include "G4VisCommandsScene.hh"

#include "G4VisManager.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4Run.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4ApplicationState.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4ios.hh"
#include <sstream>

#define G4warn G4cout

G4VVisCommandScene::G4VVisCommandScene () {}

G4VVisCommandScene::~G4VVisCommandScene () {}

G4String G4VVisCommandScene::CurrentSceneName () {
  const G4Scene* pScene = fpVisManager -> GetCurrentScene ();
  G4String currentSceneName = "none";
  if (pScene) currentSceneName = pScene -> GetName ();
  return currentSceneName;
}

////////////// /vis/scene/activateModel ////////////////////////////

G4VisCommandSceneActivateModel::G4VisCommandSceneActivateModel () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/activateModel", this);
  fpCommand -> SetGuidance
    ("Activate or de-activate model.");
  fpCommand -> SetGuidance
    ("Attempts to match search string to name of model - use unique sub-string.");
  fpCommand -> SetGuidance
    ("Use \"/vis/scene/list\" to see model names.");
  fpCommand -> SetGuidance
    ("If name == \"all\" (default), all models are activated.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("search-string", 's', omitable = true);
  parameter -> SetDefaultValue ("all");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("activate", 'b', omitable = true);
  parameter -> SetDefaultValue (true);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneActivateModel::~G4VisCommandSceneActivateModel () {
  delete fpCommand;
}

G4String G4VisCommandSceneActivateModel::GetCurrentValue(G4UIcommand*) {
  return "";
}

void G4VisCommandSceneActivateModel::SetNewValue (G4UIcommand*,
						     G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String searchString, activateString;
  std::istringstream is (newValue);
  is >> searchString >> activateString;
  G4bool activate = G4UIcommand::ConvertToBool(activateString);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4VSceneHandler* pSceneHandler = fpVisManager->GetCurrentSceneHandler();
  if (!pSceneHandler) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current sceneHandler.  Please create one." << G4endl;
    }
    return;
  }

  if (searchString == "all" && !activate) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn <<
	"WARNING: You are not allowed to de-activate all models."
	"\n  Command ignored."
	     << G4endl;
    }
    return;
  }

  G4bool any = false;

  std::vector<G4Scene::Model>& runDurationModelList =
    pScene->SetRunDurationModelList();
  for (size_t i = 0; i < runDurationModelList.size(); i++) {
    const G4String& modelName =
      runDurationModelList[i].fpModel->GetGlobalDescription();
    if (searchString == "all" || modelName.find(searchString)
	!= std::string::npos) {
      any = true;
      runDurationModelList[i].fActive = activate;
      if (verbosity >= G4VisManager::warnings) {
	G4warn << "Model \"" << modelName;
	if (activate) G4warn << "\" activated.";
	else  G4warn << "\" de-activated.";
	G4warn << G4endl;
      }
    }
  }

  std::vector<G4Scene::Model>& endOfEventModelList =
    pScene->SetEndOfEventModelList();
  for (size_t i = 0; i < endOfEventModelList.size(); i++) {
    const G4String& modelName =
      endOfEventModelList[i].fpModel->GetGlobalDescription();
    if (searchString == "all" || modelName.find(searchString)
	!= std::string::npos) {
      any = true;
      endOfEventModelList[i].fActive = activate;
      if (verbosity >= G4VisManager::warnings) {
	G4warn << "Model \"" << modelName;
	if (activate) G4warn << "\" activated.";
	else  G4warn << "\" de-activated.";
	G4warn << G4endl;
      }
    }
  }

  std::vector<G4Scene::Model>& endOfRunModelList =
    pScene->SetEndOfRunModelList();
  for (size_t i = 0; i < endOfRunModelList.size(); i++) {
    const G4String& modelName =
      endOfRunModelList[i].fpModel->GetGlobalDescription();
    if (searchString == "all" || modelName.find(searchString)
	!= std::string::npos) {
      any = true;
      endOfRunModelList[i].fActive = activate;
      if (verbosity >= G4VisManager::warnings) {
	G4warn << "Model \"" << modelName;
	if (activate) G4warn << "\" activated.";
	else  G4warn << "\" de-activated.";
	G4warn << G4endl;
      }
    }
  }

  if (!any) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "WARNING: No match found." << G4endl;
    }
    return;
  }

  CheckSceneAndNotifyHandlers (pScene);
}

////////////// /vis/scene/create ///////////////////////////////////////

G4VisCommandSceneCreate::G4VisCommandSceneCreate (): fId (0) {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString ("/vis/scene/create", this);
  fpCommand -> SetGuidance
    ("Creates an empty scene.");
  fpCommand -> SetGuidance
    ("Invents a name if not supplied.  This scene becomes current.");
  fpCommand -> SetParameterName ("scene-name", omitable = true);
}

G4VisCommandSceneCreate::~G4VisCommandSceneCreate () {
  delete fpCommand;
}

G4String G4VisCommandSceneCreate::NextName () {
  std::ostringstream oss;
  oss << "scene-" << fId;
  return oss.str();
}

G4String G4VisCommandSceneCreate::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneCreate::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String& newName = newValue;
  G4String nextName = NextName ();

  if (newName == "") {
    newName = nextName;
  }
  if (newName == nextName) fId++;

  G4SceneList& sceneList = fpVisManager -> SetSceneList ();
  std::size_t iScene, nScenes = sceneList.size ();
  for (iScene = 0; iScene < nScenes; ++iScene) {
    if (sceneList [iScene] -> GetName () == newName) break;
  }
  if (iScene < nScenes) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "WARNING: Scene \"" << newName << "\" already exists."
	     << "\n  New scene not created."
	     << G4endl;
    }
  } else {

    // Add empty scene data object to list...
    G4Scene* pScene = new G4Scene (newName);
    sceneList.push_back (pScene);
    fpVisManager -> SetCurrentScene (pScene);

    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "New empty scene \"" << newName << "\" created." << G4endl;
    }
  }
}

////////////// /vis/scene/endOfEventAction ////////////////////////////

G4VisCommandSceneEndOfEventAction::G4VisCommandSceneEndOfEventAction () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/endOfEventAction", this);
  fpCommand -> SetGuidance
    ("Accumulate or refresh the viewer for each new event.");
  fpCommand -> SetGuidance
    ("\"accumulate\": viewer accumulates hits, etc., event by event, or");
  fpCommand -> SetGuidance
    ("\"refresh\": viewer shows them at end of event or, for direct-screen"
     "\n  viewers, refreshes the screen just before drawing the next event.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("action", 's', omitable = true);
  parameter -> SetParameterCandidates ("accumulate refresh");
  parameter -> SetDefaultValue ("refresh");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("maxNumber", 'i', omitable = true);
  parameter -> SetDefaultValue (100);
  parameter -> SetGuidance
  ("Maximum number of events kept.  Unlimited if negative.");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneEndOfEventAction::~G4VisCommandSceneEndOfEventAction () {
  delete fpCommand;
}

G4String G4VisCommandSceneEndOfEventAction::GetCurrentValue(G4UIcommand*) {
  return "";
}

void G4VisCommandSceneEndOfEventAction::SetNewValue (G4UIcommand*,
						     G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String action;
  G4int maxNumberOfKeptEvents;
  std::istringstream is (newValue);
  is >> action >> maxNumberOfKeptEvents;

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4VSceneHandler* pSceneHandler = fpVisManager->GetCurrentSceneHandler();
  if (!pSceneHandler) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current sceneHandler.  Please create one." << G4endl;
    }
    return;
  }

  if (action == "accumulate") {
    pScene->SetRefreshAtEndOfEvent(false);
    pScene->SetMaxNumberOfKeptEvents(maxNumberOfKeptEvents);
  }
  else if (action == "refresh") {
    if (!pScene->GetRefreshAtEndOfRun()) {
      if (verbosity >= G4VisManager::errors) {
	G4warn <<
	  "ERROR: Cannot refresh events unless runs refresh too."
	  "\n  Use \"/vis/scene/endOfRun refresh\"."
	       << G4endl;
      }
    } else {
      pScene->SetRefreshAtEndOfEvent(true);
      pScene->SetMaxNumberOfKeptEvents(maxNumberOfKeptEvents);
      pSceneHandler->SetMarkForClearingTransientStore(true);
    }
  }
  else {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
	"ERROR: unrecognised parameter \"" << action << "\"."
             << G4endl;
    }
    return;
  }

  // Change of transients behaviour, so...
  fpVisManager->ResetTransientsDrawnFlags();

  // Are there any events currently kept...
  size_t nCurrentlyKept    = 0;
  G4RunManager* runManager = G4RunManagerFactory::GetMasterRunManager();
  if(runManager)
  {
    const G4Run* currentRun = runManager->GetCurrentRun();
    if(currentRun)
    {
      const std::vector<const G4Event*>* events = currentRun->GetEventVector();
      if(events)
        nCurrentlyKept = events->size();
    }
  }

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "End of event action set to ";
    if (pScene->GetRefreshAtEndOfEvent()) G4cout << "\"refresh\".";
    else {
      G4cout << "\"accumulate\"."
	"\n  Maximum number of events to be kept: "
	     << maxNumberOfKeptEvents
	     << " (unlimited if negative)."
	"\n  This may be changed with, e.g., "
	"\"/vis/scene/endOfEventAction accumulate 1000\".";
    }
    G4cout << G4endl;
  }

  if (!pScene->GetRefreshAtEndOfEvent() &&
      maxNumberOfKeptEvents != 0 &&
      verbosity >= G4VisManager::warnings) {
    G4warn << "WARNING: ";
    if (nCurrentlyKept) {
      G4warn <<
	"\n  There are currently " << nCurrentlyKept
	     << " events kept for refreshing and/or reviewing.";
    } else {
      G4warn << "The vis manager will keep ";
      if (maxNumberOfKeptEvents < 0) G4warn << "an unlimited number of";
      else G4warn << "up to " << maxNumberOfKeptEvents;
      G4warn << " events.";
      if (maxNumberOfKeptEvents > 1 || maxNumberOfKeptEvents < 0)
	G4warn <<
	  "\n  This may use a lot of memory."
	  "\n  It may be changed with, e.g., "
	  "\"/vis/scene/endOfEventAction accumulate 10\".";
    }
    G4warn << G4endl;
  }
}

////////////// /vis/scene/endOfRunAction ////////////////////////////

G4VisCommandSceneEndOfRunAction::G4VisCommandSceneEndOfRunAction () {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString ("/vis/scene/endOfRunAction", this);
  fpCommand -> SetGuidance
    ("Accumulate or refresh the viewer for each new run.");
  fpCommand -> SetGuidance
    ("\"accumulate\": viewer accumulates hits, etc., run by run, or");
  fpCommand -> SetGuidance
    ("\"refresh\": viewer shows them at end of run or, for direct-screen"
     "\n  viewers, refreshes the screen just before drawing the first"
     "\n  event of the next run.");
  fpCommand -> SetGuidance ("The detector remains or is redrawn.");
  fpCommand -> SetParameterName ("action", omitable = true);
  fpCommand -> SetCandidates ("accumulate refresh");
  fpCommand -> SetDefaultValue ("refresh");
}

G4VisCommandSceneEndOfRunAction::~G4VisCommandSceneEndOfRunAction () {
  delete fpCommand;
}

G4String G4VisCommandSceneEndOfRunAction::GetCurrentValue(G4UIcommand*) {
  return "";
}

void G4VisCommandSceneEndOfRunAction::SetNewValue (G4UIcommand*,
						     G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String action;
  std::istringstream is (newValue);
  is >> action;

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4VSceneHandler* pSceneHandler = fpVisManager->GetCurrentSceneHandler();
  if (!pSceneHandler) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current sceneHandler.  Please create one." << G4endl;
    }
    return;
  }

  if (action == "accumulate") {
    if (pScene->GetRefreshAtEndOfEvent()) {
      if (verbosity >= G4VisManager::errors) {
	G4warn <<
	  "ERROR: Cannot accumulate runs unless events accumulate too."
	  "\n  Use \"/vis/scene/endOfEventAction accumulate\"."
	       << G4endl;
      }
    }
    else {
      pScene->SetRefreshAtEndOfRun(false);
    }
  }
  else if (action == "refresh") {
    pScene->SetRefreshAtEndOfRun(true);
    pSceneHandler->SetMarkForClearingTransientStore(true);
  }
  else {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
	"ERROR: unrecognised parameter \"" << action << "\"."
             << G4endl;
    }
    return;
  }

  // Change of transients behaviour, so...
  fpVisManager->ResetTransientsDrawnFlags();

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "End of run action set to \"";
    if (pScene->GetRefreshAtEndOfRun()) G4cout << "refresh";
    else G4cout << "accumulate";
    G4cout << "\"" << G4endl;
  }
}

////////////// /vis/scene/list ///////////////////////////////////////

G4VisCommandSceneList::G4VisCommandSceneList () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/list", this);
  fpCommand -> SetGuidance ("Lists scene(s).");
  fpCommand -> SetGuidance
    ("\"help /vis/verbose\" for definition of verbosity.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("scene-name", 's', omitable = true);
  parameter -> SetDefaultValue ("all");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("verbosity", 's', omitable = true);
  parameter -> SetDefaultValue ("warnings");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneList::~G4VisCommandSceneList () {
  delete fpCommand;
}

G4String G4VisCommandSceneList::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneList::SetNewValue (G4UIcommand*, G4String newValue) {
  G4String name, verbosityString;
  std::istringstream is (newValue);
  is >> name >> verbosityString;
  G4VisManager::Verbosity verbosity =
    fpVisManager->GetVerbosityValue(verbosityString);
  const G4Scene* currentScene = fpVisManager -> GetCurrentScene ();
  G4String currentName;
  if (currentScene) currentName = currentScene->GetName();

  G4SceneList& sceneList = fpVisManager -> SetSceneList ();
  std::size_t iScene, nScenes = sceneList.size ();
  G4bool found = false;
  for (iScene = 0; iScene < nScenes; ++iScene) {
    G4Scene* pScene = sceneList [iScene];
    const G4String& iName = pScene -> GetName ();
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
    G4cout << " scene \"" << iName << "\"";
    if (verbosity >= G4VisManager::warnings) {
      std::size_t i;
      G4cout << "\n  Run-duration models:";
      std::size_t nRunModels = pScene -> GetRunDurationModelList ().size ();
      if (nRunModels == 0) {
	G4cout << " none.";
      }
      for (i = 0; i < nRunModels; ++i) {
	if (pScene -> GetRunDurationModelList()[i].fActive)
	  G4cout << "\n   Active:   ";
	else G4cout << "\n   Inactive: ";
	G4VModel* pModel = pScene -> GetRunDurationModelList()[i].fpModel;
	G4cout << pModel -> GetGlobalDescription ();
      }
      G4cout << "\n  End-of-event models:";
      std::size_t nEOEModels = pScene -> GetEndOfEventModelList ().size ();
      if (nEOEModels == 0) {
	G4cout << " none.";
      }
      for (i = 0; i < nEOEModels; ++i) {
	if (pScene -> GetEndOfEventModelList()[i].fActive)
	  G4cout << "\n   Active:   ";
	else G4cout << "\n   Inactive: ";
	G4VModel* pModel = pScene -> GetEndOfEventModelList()[i].fpModel;
	G4cout << pModel -> GetGlobalDescription ();
      }
      G4cout << "\n  End-of-run models:";
      std::size_t nEORModels = pScene -> GetEndOfRunModelList ().size ();
      if (nEORModels == 0) {
	G4cout << " none.";
      }
      for (i = 0; i < nEORModels; ++i) {
	if (pScene -> GetEndOfRunModelList()[i].fActive)
	  G4cout << "\n   Active:   ";
	else G4cout << "\n   Inactive: ";
	G4VModel* pModel = pScene -> GetEndOfRunModelList()[i].fpModel;
	G4cout << pModel -> GetGlobalDescription ();
      }
    }
    if (verbosity >= G4VisManager::parameters) {
      G4cout << "\n  " << *sceneList [iScene];
    }
    G4cout << G4endl;
  }
  if (!found) {
    G4warn << "No scenes found";
    if (name != "all") {
      G4warn << " of name \"" << name << "\"";
    }
    G4warn << "." << G4endl;
  }
}

////////////// /vis/scene/notifyHandlers /////////////////////////

G4VisCommandSceneNotifyHandlers::G4VisCommandSceneNotifyHandlers () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/notifyHandlers", this);
  fpCommand -> SetGuidance
    ("Notifies scene handlers and forces re-rendering.");
  fpCommand -> SetGuidance
    ("Notifies the handler(s) of the specified scene and forces a"
     "\nreconstruction of any graphical databases."
     "\nClears and refreshes all viewers of current scene."
     "\n  The default action \"refresh\" does not issue \"update\" (see"
     "\n    /vis/viewer/update)."
     "\nIf \"flush\" is specified, it issues an \"update\" as well as"
     "\n  \"refresh\" - \"update\" and initiates post-processing"
     "\n  for graphics systems which need it.");
  fpCommand -> SetGuidance 
    ("The default for <scene-name> is the current scene name.");
  fpCommand -> SetGuidance
    ("This command does not change current scene, scene handler or viewer.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("scene-name", 's',
				 omitable = true);
  parameter -> SetCurrentAsDefault(true);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("refresh-flush", 's',
				 omitable = true);
  parameter -> SetDefaultValue("refresh");
  parameter -> SetParameterCandidates("r refresh f flush");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneNotifyHandlers::~G4VisCommandSceneNotifyHandlers () {
  delete fpCommand;
}

G4String G4VisCommandSceneNotifyHandlers::GetCurrentValue(G4UIcommand*) {
  return CurrentSceneName ();
}

void G4VisCommandSceneNotifyHandlers::SetNewValue (G4UIcommand*,
						   G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String sceneName, refresh_flush;
  std::istringstream is (newValue);
  is >> sceneName >> refresh_flush;
  G4bool flush = false;
  if (refresh_flush[0] == 'f') flush = true;

  const G4SceneList& sceneList = fpVisManager -> GetSceneList ();
  G4SceneHandlerList& sceneHandlerList =
    fpVisManager -> SetAvailableSceneHandlers ();

  // Check scene name.
  const std::size_t nScenes = sceneList.size ();
  std::size_t iScene;
  for (iScene = 0; iScene < nScenes; ++iScene) {
    G4Scene* scene = sceneList [iScene];
    if (sceneName == scene -> GetName ()) break;
  }
  if (iScene >= nScenes ) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "WARNING: Scene \"" << sceneName << "\" not found."
	"\n  /vis/scene/list to see scenes."
	     << G4endl;
    }
    return;
  }

  // Store current context...
  G4VSceneHandler* pCurrentSceneHandler =
    fpVisManager -> GetCurrentSceneHandler();
  if (!pCurrentSceneHandler) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "WARNING: No current scene handler."
	     << G4endl;
    }
    return;
  }
  G4VViewer* pCurrentViewer = fpVisManager -> GetCurrentViewer();
  if (!pCurrentViewer) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "WARNING: No current viewer."
	     << G4endl;
    }
    return;
  }
  G4Scene* pCurrentScene = fpVisManager -> GetCurrentScene();
  if (!pCurrentScene) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "WARNING: No current scene."
	     << G4endl;
    }
    return;
  }

  G4VisManager::Verbosity currentVerbosity = fpVisManager -> GetVerbosity();

  // Suppress messages during this process (only print errors)...
  //fpVisManager -> SetVerboseLevel(G4VisManager::errors);

  // For each scene handler, if it contains the scene, clear and
  // rebuild the graphical database, then for each viewer set (make
  // current), clear, (re)draw, and show.
  const std::size_t nSceneHandlers = sceneHandlerList.size ();
  for (std::size_t iSH = 0; iSH < nSceneHandlers; ++iSH) {
    G4VSceneHandler* aSceneHandler = sceneHandlerList [iSH];
    G4Scene* aScene = aSceneHandler -> GetScene ();
    if (aScene) {
      const G4String& aSceneName = aScene -> GetName ();
      if (sceneName == aSceneName) {
	aScene->CalculateExtent();  // Check and recalculate extent
	G4ViewerList& viewerList = aSceneHandler -> SetViewerList ();
	const std::size_t nViewers = viewerList.size ();
	for (std::size_t iV = 0; iV < nViewers; ++iV) {
	  G4VViewer* aViewer = viewerList [iV];
	  // Force rebuild of graphical database, if any.
	  aViewer -> NeedKernelVisit();
	  if (aViewer->GetViewParameters().IsAutoRefresh()) {
	    aSceneHandler -> SetCurrentViewer (aViewer);
	    // Ensure consistency of vis manager...
	    fpVisManager -> SetCurrentViewer(aViewer);
	    fpVisManager -> SetCurrentSceneHandler(aSceneHandler);
	    fpVisManager -> SetCurrentScene(aScene);
	    aViewer -> SetView ();
	    aViewer -> ClearView ();
	    aViewer -> DrawView ();
	    if (flush) aViewer -> ShowView ();
	    if (verbosity >= G4VisManager::confirmations) {
	      G4cout << "Viewer \"" << aViewer -> GetName ()
		     << "\" of scene handler \"" << aSceneHandler -> GetName ()
		     << "\"\n  ";
	      if (flush) G4cout << "flushed";
	      else G4cout << "refreshed";
	      G4cout << " at request of scene \"" << sceneName
		     << "\"." << G4endl;
	    }
	  } else {
	    if (verbosity >= G4VisManager::confirmations) {
	      G4cout << "NOTE: The scene, \""
		     << sceneName
		     << "\", of viewer \""
		     << aViewer -> GetName ()
		     << "\"\n  of scene handler \""
		     << aSceneHandler -> GetName ()
		     << "\"  has changed.  To see effect,"
		     << "\n  \"/vis/viewer/select "
		     << aViewer -> GetShortName ()
		     << "\" and \"/vis/viewer/rebuild\"."
		     << G4endl;
	    }
	  }
	}
      }
    }
    else {
      if (verbosity >= G4VisManager::warnings) {
	G4warn << "WARNING: G4VisCommandSceneNotifyHandlers: scene handler \""
	       << aSceneHandler->GetName()
	       << "\" has a null scene."
	       << G4endl;
      }
    }
  }

  // Reclaim original context - but set viewer first, then scene
  // handler, because the latter might have been created very recently
  // and, not yet having a viewer, the current viewer will,
  // temporarily, refer to another scene handler.  SetCurrentViewer
  // actually resets the scene handler, which is what we don't want,
  // so we set it again on the next line...
  fpVisManager -> SetCurrentViewer(pCurrentViewer);
  fpVisManager -> SetCurrentSceneHandler(pCurrentSceneHandler);
  fpVisManager -> SetCurrentScene(pCurrentScene);
  fpVisManager -> SetVerboseLevel(currentVerbosity);
  // Take care of special case of scene handler with no viewer yet.  
  if (pCurrentSceneHandler) {
    G4ViewerList& viewerList = pCurrentSceneHandler -> SetViewerList ();
    const std::size_t nViewers = viewerList.size ();
    if (nViewers) {
      pCurrentSceneHandler -> SetCurrentViewer (pCurrentViewer);
      // JA: I don't think we need this. SetView will be called when needed.
      // if (pCurrentViewer && pCurrentSceneHandler->GetScene()) {
      //   pCurrentViewer -> SetView ();
      // }
    }
  }
}

////////////// /vis/scene/removeModel ////////////////////////////

G4VisCommandSceneRemoveModel::G4VisCommandSceneRemoveModel () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/removeModel", this);
  fpCommand -> SetGuidance("Remove model.");
  fpCommand -> SetGuidance
  ("Attempts to match search string to name of model - use unique sub-string.");
  fpCommand -> SetGuidance
  ("Use \"/vis/scene/list\" to see model names.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("search-string", 's', omitable = false);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneRemoveModel::~G4VisCommandSceneRemoveModel () {
  delete fpCommand;
}

G4String G4VisCommandSceneRemoveModel::GetCurrentValue(G4UIcommand*) {
  return "";
}

void G4VisCommandSceneRemoveModel::SetNewValue (G4UIcommand*,
						  G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String searchString;
  std::istringstream is (newValue);
  is >> searchString;

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4VSceneHandler* pSceneHandler = fpVisManager->GetCurrentSceneHandler();
  if (!pSceneHandler) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current sceneHandler.  Please create one." << G4endl;
    }
    return;
  }

  G4bool any = false;

  std::vector<G4Scene::Model>& runDurationModelList =
  pScene->SetRunDurationModelList();
  for (size_t i = 0; i < runDurationModelList.size(); i++) {
    const G4String& modelName =
    runDurationModelList[i].fpModel->GetGlobalDescription();
    if (modelName.find(searchString) != std::string::npos) {
      runDurationModelList.erase(runDurationModelList.begin()+i);
      any = true;
      if (verbosity >= G4VisManager::warnings) {
	G4warn << "Model \"" << modelName << "\" removed." << G4endl;
      }
      break;  // Allow only one model at a time to be removed.
    }
  }

  std::vector<G4Scene::Model>& endOfEventModelList =
  pScene->SetEndOfEventModelList();
  for (size_t i = 0; i < endOfEventModelList.size(); i++) {
    const G4String& modelName =
    endOfEventModelList[i].fpModel->GetGlobalDescription();
    if (modelName.find(searchString) != std::string::npos) {
      endOfEventModelList.erase(endOfEventModelList.begin()+i);
      any = true;
      if (verbosity >= G4VisManager::warnings) {
	G4warn << "Model \"" << modelName << "\" removed." << G4endl;
      }
      break;  // Allow only one model at a time to be removed.
    }
  }

  std::vector<G4Scene::Model>& endOfRunModelList =
  pScene->SetEndOfRunModelList();
  for (size_t i = 0; i < endOfRunModelList.size(); i++) {
    const G4String& modelName =
    endOfRunModelList[i].fpModel->GetGlobalDescription();
    if (modelName.find(searchString) != std::string::npos) {
      endOfRunModelList.erase(endOfRunModelList.begin()+i);
      any = true;
      if (verbosity >= G4VisManager::warnings) {
	G4warn << "Model \"" << modelName << "\" removed." << G4endl;
      }
      break;  // Allow only one model at a time to be removed.
    }
  }

  if (!any) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "WARNING: No match found." << G4endl;
    }
    return;
  }

  CheckSceneAndNotifyHandlers (pScene);
}

////////////// /vis/scene/select ///////////////////////////////////////

G4VisCommandSceneSelect::G4VisCommandSceneSelect () {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString ("/vis/scene/select", this);
  fpCommand -> SetGuidance ("Selects a scene");
  fpCommand -> SetGuidance
  ("Makes the scene current.  \"/vis/scene/list\" to see"
   "\n possible scene names.");
  fpCommand -> SetParameterName ("scene-name", omitable = false);
}

G4VisCommandSceneSelect::~G4VisCommandSceneSelect () {
  delete fpCommand;
}

G4String G4VisCommandSceneSelect::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneSelect::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String& selectName = newValue;
  G4SceneList& sceneList = fpVisManager -> SetSceneList ();
  std::size_t iScene, nScenes = sceneList.size ();
  for (iScene = 0; iScene < nScenes; ++iScene) {
    if (sceneList [iScene] -> GetName () == selectName) break;
  }
  if (iScene >= nScenes) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "WARNING: Scene \"" << selectName
      << "\" not found - \"/vis/scene/list\" to see possibilities."
      << G4endl;
    }
    return;
  }

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Scene \"" << selectName
    << "\" selected." << G4endl;
  }

  CheckSceneAndNotifyHandlers (sceneList [iScene]);
}

////////////// /vis/scene/showExtents ///////////////////////////////////////

G4VisCommandSceneShowExtents::G4VisCommandSceneShowExtents () {
  fpCommand = new G4UIcmdWithoutParameter ("/vis/scene/showExtents", this);
  fpCommand -> SetGuidance ("Prints and draws extents of models in a scene");
}

G4VisCommandSceneShowExtents::~G4VisCommandSceneShowExtents () {
  delete fpCommand;
}

G4String G4VisCommandSceneShowExtents::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneShowExtents::SetNewValue (G4UIcommand*, G4String) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VSceneHandler* pCurrentSceneHandler =
  fpVisManager -> GetCurrentSceneHandler();
  if (!pCurrentSceneHandler) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "WARNING: No current scene handler."
      << G4endl;
    }
    return;
  }
  G4VViewer* pCurrentViewer = fpVisManager -> GetCurrentViewer();
  if (!pCurrentViewer) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "WARNING: No current viewer."
      << G4endl;
    }
    return;
  }
  G4Scene* pCurrentScene = fpVisManager -> GetCurrentScene();
  if (!pCurrentScene) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn << "WARNING: No current scene."
      << G4endl;
    }
    return;
  }

  G4cout << "\n  Run-duration models:";
  std::size_t nRunModels = pCurrentScene -> GetRunDurationModelList ().size ();
  if (nRunModels == 0) {
    G4cout << " none.";
  }
  for (std::size_t i = 0; i < nRunModels; ++i) {
    if (pCurrentScene -> GetRunDurationModelList()[i].fActive)
      G4cout << "\n   Active:   ";
    else G4cout << "\n   Inactive: ";
    G4VModel* pModel = pCurrentScene -> GetRunDurationModelList()[i].fpModel;
    const G4VisExtent& transformedExtent = pModel -> GetExtent();
    G4cout << pModel -> GetGlobalDescription ()
    << "\n" << transformedExtent;
    DrawExtent(transformedExtent);
  }
  G4cout << "\n  End-of-event models:";
  std::size_t nEOEModels = pCurrentScene -> GetEndOfEventModelList ().size ();
  if (nEOEModels == 0) {
    G4cout << " none.";
  }
  for (std::size_t i = 0; i < nEOEModels; ++i) {
    if (pCurrentScene -> GetEndOfEventModelList()[i].fActive)
      G4cout << "\n   Active:   ";
    else G4cout << "\n   Inactive: ";
    G4VModel* pModel = pCurrentScene -> GetEndOfEventModelList()[i].fpModel;
    const G4VisExtent& transformedExtent = pModel -> GetExtent();
    G4cout << pModel -> GetGlobalDescription ()
    << "\n" << transformedExtent;
    DrawExtent(transformedExtent);
  }
  G4cout << "\n  End-of-run models:";
  std::size_t nEORModels = pCurrentScene -> GetEndOfRunModelList ().size ();
  if (nEORModels == 0) {
    G4cout << " none.";
  }
  for (std::size_t i = 0; i < nEORModels; ++i) {
    if (pCurrentScene -> GetEndOfRunModelList()[i].fActive)
      G4cout << "\n   Active:   ";
    else G4cout << "\n   Inactive: ";
    G4VModel* pModel = pCurrentScene -> GetEndOfRunModelList()[i].fpModel;
    const G4VisExtent& transformedExtent = pModel -> GetExtent();
    G4cout << pModel -> GetGlobalDescription ()
    << "\n" << transformedExtent;
    DrawExtent(transformedExtent);
  }
  G4cout << "\n  Overall extent:\n";
  DrawExtent(pCurrentScene->GetExtent());
  G4cout << G4endl;
}
