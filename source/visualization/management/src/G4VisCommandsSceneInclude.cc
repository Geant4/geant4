// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsSceneInclude.cc,v 1.5 2000-05-19 09:18:01 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/scene commands - John Allison  9th August 1998

#include "G4VisCommandsSceneInclude.hh"

#include "G4VisManager.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4HitsModel.hh"
#include "G4TrajectoriesModel.hh"

////////////// /vis/scene/include/hits ///////////////////////////////////////

G4VisCommandSceneIncludeHits::G4VisCommandSceneIncludeHits () {
  fpCommand = new G4UIcmdWithoutParameter ("/vis/scene/include/hits", this);
  fpCommand -> AvailableForStates (Idle, GeomClosed);
  fpCommand -> SetGuidance
    ("Deprecated command - please use /vis/scene/add/hits.");
  fpCommand -> SetGuidance
    ("Includes hits in current scene.");
  fpCommand -> SetGuidance
    ("Hits are drawn at end of event when the scene in which"
     " they are included is current.");
}

G4VisCommandSceneIncludeHits::~G4VisCommandSceneIncludeHits () {
  delete fpCommand;
}

G4String G4VisCommandSceneIncludeHits::GetCurrentValue (G4UIcommand* command) {
  return "";
}

void G4VisCommandSceneIncludeHits::SetNewValue (G4UIcommand* command,
						G4String newValue) {
  G4cout <<
    "/vis/scene/include/hits is a deprecated command - please use"
    " /vis/scene/add/hits."
	 << G4endl;

  G4SceneList& list = fpVisManager -> SetSceneList ();
  if (list.isEmpty ()) {
    G4cout << "No scenes - please create one before adding anything."
	   << G4endl;
    return;
  }

  G4HitsModel* model = new G4HitsModel;
  G4Scene* pCurrentScene = fpVisManager -> GetCurrentScene ();
  const G4String& currentSceneName = pCurrentScene -> GetName ();
  pCurrentScene -> AddEndOfEventModel (model);
  G4cout << "Hits will be drawn in scene \""
	 << currentSceneName << "\"."
	 << G4endl;
}

////////////// /vis/scene/include/trajectories ///////////////////////////////////////

G4VisCommandSceneIncludeTrajectories::G4VisCommandSceneIncludeTrajectories () {
  fpCommand = new G4UIcmdWithoutParameter
    ("/vis/scene/include/trajectories", this);
  fpCommand -> SetGuidance
    ("Deprecated command - please use /vis/scene/add/trajectories.");
  fpCommand -> AvailableForStates (Idle, GeomClosed);
  fpCommand -> SetGuidance
    ("Includes trajectories in current scene.");
  fpCommand -> SetGuidance
    ("Trajectories are drawn at end of event when the scene in which"
     " they are included is current.");
}

G4VisCommandSceneIncludeTrajectories::~G4VisCommandSceneIncludeTrajectories () {
  delete fpCommand;
}

G4String G4VisCommandSceneIncludeTrajectories::GetCurrentValue (G4UIcommand* command) {
  return "";
}

void G4VisCommandSceneIncludeTrajectories::SetNewValue (G4UIcommand* command,
					      G4String newValue) {
  G4cout <<
    "/vis/scene/include/trajectories is a deprecated command - please use"
    "\n  /vis/scene/add/trajectories."
	 << G4endl;

  G4SceneList& list = fpVisManager -> SetSceneList ();
  if (list.isEmpty ()) {
    G4cout << "No scenes - please create one before adding anything."
	   << G4endl;
    return;
  }

  G4TrajectoriesModel* model = new G4TrajectoriesModel;
  G4Scene* pCurrentScene = fpVisManager -> GetCurrentScene ();
  const G4String& currentSceneName = pCurrentScene -> GetName ();
  pCurrentScene -> AddEndOfEventModel (model);
  G4cout << "Trajectories will be drawn in scene \""
	 << currentSceneName << "\"."
	 << G4endl;
}
