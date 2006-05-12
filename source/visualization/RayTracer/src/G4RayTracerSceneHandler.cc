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
// $Id: G4RayTracerSceneHandler.cc,v 1.7 2006-05-12 12:43:43 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4RayTracerSceneHandler.hh"
#include "G4VisManager.hh"

G4RayTracerSceneHandler::G4RayTracerSceneHandler(G4VGraphicsSystem& system,
						 const G4String& name):
  G4VSceneHandler(system, fSceneIdCount++, name)
{
  // Ray Tracer does not use the vis concept of scene (G4Scene) - keep
  // vis manager happy when someone opens a ray tracer with "/vis/open
  // RayTracer" but uses the ray tracer with "/vis/rayTracer" commands
  // before creating any scenes, for example, instead of using
  // "/vis/drawVolume"...
  G4VisManager* visManager = G4VisManager::GetInstance();
  if(visManager) {
    G4Scene* pScene = visManager->GetCurrentScene();
    if (!pScene) {
      // Create new scene like /vis/scene/create...
      fpScene = new G4Scene("dummy-ray-tracer-scene");
      // Avoid code triggered at end of events...
      fpScene->SetRefreshAtEndOfEvent(false);
      // Avoid re-computing transients.
      fpScene->SetRecomputeTransients(false);
      // Add to vis manager list; ownership thereby passes to vis manager...
      visManager->SetSceneList().push_back(fpScene);
      // ...and make current...
      visManager->SetCurrentScene(fpScene);
    }
  }
}

G4RayTracerSceneHandler::~G4RayTracerSceneHandler()
{}

G4int G4RayTracerSceneHandler::fSceneIdCount = 0;
