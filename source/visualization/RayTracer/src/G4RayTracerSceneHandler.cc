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
// $Id: G4RayTracerSceneHandler.cc,v 1.4 2006/01/11 18:01:33 allison Exp $
// GEANT4 tag $Name: geant4-08-00-patch-01 $

#include "G4RayTracerSceneHandler.hh"
#include "G4VisManager.hh"
#include "G4NullModel.hh"

G4Scene G4RayTracerSceneHandler::fDummyRayTracerScene
("G4RayTracerSceneHandler::fDummyRayTracerScene");

G4RayTracerSceneHandler::G4RayTracerSceneHandler(G4VGraphicsSystem& system,
						 const G4String& name):
  G4VSceneHandler(system, fSceneIdCount++, name)
{
  // Ray Tracer does not use the vis concept of scene (G4Scene) - keep
  // vis manager happy when someone opens a ray tracer with "/vis/open
  // RayTracer" but uses the ray tracer with "/vis/rayTracer" commands
  // before creating any scenes, for example, instead of using
  // "/vis/drawVolume"...
  fpScene = &fDummyRayTracerScene;
  if (fpScene->IsEmpty())
    fpScene->SetRunDurationModelList().push_back(new G4NullModel);
  G4VisManager* visManager = G4VisManager::GetInstance();
  if(visManager) {
    G4Scene* pScene = visManager->GetCurrentScene();
    if (!pScene) {
      visManager->SetCurrentScene(fpScene);
    }
  }
}

G4RayTracerSceneHandler::~G4RayTracerSceneHandler()
{}

G4int G4RayTracerSceneHandler::fSceneIdCount = 0;
