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
// $Id: G4RayTracerSceneHandler.cc,v 1.8 2006-06-29 21:24:13 gunter Exp $
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
