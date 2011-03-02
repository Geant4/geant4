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
// $Id: G4RayTracerXViewer.cc,v 1.8 2009-09-16 16:56:52 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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
{}

G4RayTracerXViewer::~G4RayTracerXViewer() {}

void G4RayTracerXViewer::Initialise() {

  if (!theTracer) {
    G4cerr << "G4RayTracerXViewer::Initialise: No tracer" << G4endl;
    fViewId = -1;  // This flags an error.
    return;
  }

  theTracer->SetNColumn(fVP.GetWindowSizeHintX());
  theTracer->SetNRow(fVP.GetWindowSizeHintY());

  // Set up X Window...
  G4RTXScanner* theXScanner = (G4RTXScanner*)theTracer->GetScanner();
  if (!theXScanner->GetXWindow(fName,fVP)) {
    G4cerr << "G4RayTracerXViewer::Initialise: No scanner" << G4endl;
    fViewId = -1;  // This flags an error.
    return;
  }
}

#endif
