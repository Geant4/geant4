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
// $Id: G4RayTracer.cc 103626 2017-04-19 13:29:18Z gcosmo $

#include "G4RayTracer.hh"
#include "G4RayTracerFeatures.hh"
#include "G4RayTracerSceneHandler.hh"
#include "G4RayTracerViewer.hh"
#ifdef G4MULTITHREADED
#include "G4TheMTRayTracer.hh"
#else
#include "G4TheRayTracer.hh"
#endif

G4RayTracer::G4RayTracer():
  G4VGraphicsSystem("RayTracer",
		     "RayTracer",
		     RAYTRACER_FEATURES,
		     G4VGraphicsSystem::threeD)
{
#ifdef G4MULTITHREADED
  theRayTracer = new G4TheMTRayTracer;  // Establish default ray tracer.
#else
  theRayTracer = new G4TheRayTracer;  // Establish default ray tracer.
#endif
}

G4RayTracer::~G4RayTracer()
{
  delete theRayTracer;
}

G4VSceneHandler* G4RayTracer::CreateSceneHandler (const G4String& name) {
  G4VSceneHandler* pScene = new G4RayTracerSceneHandler (*this, name);
  return pScene;
}

G4VViewer* G4RayTracer::CreateViewer (G4VSceneHandler& sceneHandler,
				      const G4String& name) {
  G4VViewer* pViewer = new G4RayTracerViewer
  (sceneHandler, name, theRayTracer);
  if (pViewer) {
    if (pViewer->GetViewId() < 0) {
      G4cout <<
        "G4RayTracer::CreateViewer: ERROR flagged by negative"
        " view id in G4RayTracerViewer creation."
        "\n Destroying view and returning null pointer."
             << G4endl;
      delete pViewer;
      pViewer = 0;
    }
  }
  else {
    G4cout <<
      "G4RayTracer::CreateViewer: ERROR: null pointer on new G4RayTracerViewer."
           << G4endl;
  }
  return pViewer;
}
