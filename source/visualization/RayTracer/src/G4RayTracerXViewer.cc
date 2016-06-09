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
// $Id: G4RayTracerXViewer.cc,v 1.4 2006/01/11 18:01:33 allison Exp $
// GEANT4 tag $Name: geant4-08-00-patch-01 $

#ifdef G4VIS_BUILD_RAYTRACERX_DRIVER

#include "G4RayTracerXViewer.hh"

#include "G4VSceneHandler.hh"
#include "G4Scene.hh"
#include "G4TheRayTracer.hh"
#include "G4RTJpegMaker.hh"
#include "G4RTXScanner.hh"
#include "G4UImanager.hh"

#include <cstdlib>

G4RayTracerXViewer::G4RayTracerXViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  G4RayTracerViewer(sceneHandler,
		    name,
		    new G4TheRayTracer(new G4RTJpegMaker, new G4RTXScanner))
{
  // Set up X Window...
  G4RTXScanner* theXScanner = (G4RTXScanner*)theTracer->GetScanner();
  if (!theXScanner->GetXWindow(fName,fVP)) {
    fViewId = -1;  // This flags an error.
    return;
  }
}

G4RayTracerXViewer::~G4RayTracerXViewer() {}

#endif
