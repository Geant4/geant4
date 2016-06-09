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
// $Id: G4RayTracer.cc,v 1.20 2006/01/11 18:01:33 allison Exp $
// GEANT4 tag $Name: geant4-08-00-patch-01 $

#include "G4RayTracer.hh"
#include "G4RayTracerFeatures.hh"
#include "G4RayTracerSceneHandler.hh"
#include "G4RayTracerViewer.hh"
#include "G4TheRayTracer.hh"

G4RayTracer::G4RayTracer():
  G4VGraphicsSystem("RayTracer",
		     "RayTracer",
		     RAYTRACER_FEATURES,
		     G4VGraphicsSystem::threeD)
{
  new G4TheRayTracer;  // Establish default ray tracer.
}

G4RayTracer::~G4RayTracer()
{}

G4VSceneHandler* G4RayTracer::CreateSceneHandler (const G4String& name) {
  G4VSceneHandler* pScene = new G4RayTracerSceneHandler (*this, name);
  return pScene;
}

G4VViewer* G4RayTracer::CreateViewer (G4VSceneHandler& sceneHandler,
				      const G4String& name) {
  G4VViewer* pViewer = new G4RayTracerViewer (sceneHandler, name);
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
