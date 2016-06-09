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
// $Id: G4VVisCommand.cc,v 1.17 2005/11/22 17:12:43 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $

// Base class for visualization commands - John Allison  9th August 1998
// It is really a messenger - we have one command per messenger.

#include "G4VVisCommand.hh"

#include "G4UIcommand.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include <sstream>

G4VVisCommand::~G4VVisCommand () {}

G4VisManager* G4VVisCommand::fpVisManager = 0;

G4String G4VVisCommand::ConvertToString
(G4double x, G4double y, const char * unitName)
{
  G4double uv = G4UIcommand::ValueOf(unitName);
  
  std::ostringstream oss;
  oss << x/uv << " " << y/uv << " " << unitName;
  return oss.str();
}

void G4VVisCommand::ConvertToDoublePair(const G4String& paramString,
					G4double& xval,
					G4double& yval)
{
  G4double x, y;
  char unts[30];
  
  std::istringstream is(paramString);
  is >> x >> y >> unts;
  G4String unt = unts;

  xval = x*G4UIcommand::ValueOf(unt);
  yval = y*G4UIcommand::ValueOf(unt);

  return;
}

void G4VVisCommand::UpdateVisManagerScene
(const G4String& sceneName) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  const G4SceneList& sceneList = fpVisManager -> GetSceneList ();

  G4int iScene, nScenes = sceneList.size ();
  for (iScene = 0; iScene < nScenes; iScene++) {
    if (sceneList [iScene] -> GetName () == sceneName) break;
  }

  G4Scene* pScene = 0;  // Zero unless scene has been found...
  if (iScene < nScenes) {
    pScene = sceneList [iScene];
  }

  if (!pScene) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout << "WARNING: Scene \"" << sceneName << "\" not found."
	     << G4endl;
    }
    return;
  }

  fpVisManager -> SetCurrentScene (pScene);

  // Scene has changed.  Refresh viewers of all sceneHandlers using
  // this scene...
  G4VViewer* pViewer = fpVisManager -> GetCurrentViewer();
  G4VSceneHandler* sceneHandler = fpVisManager -> GetCurrentSceneHandler();
  if (sceneHandler && sceneHandler -> GetScene ()) {
    if (pViewer) {
      G4UImanager::GetUIpointer () ->
	ApplyCommand ("/vis/scene/notifyHandlers");
    }
  }
}
