// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayTracerViewer.cc,v 1.1 2000-02-23 16:03:53 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4RayTracerViewer.hh"

#include "G4ios.hh"
#include "g4std/strstream"

#include "G4VSceneHandler.hh"
#include "G4RayTracer.hh"
#include "G4UImanager.hh"

G4RayTracerViewer::G4RayTracerViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  fFileCount(0),
  G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name) {}

G4RayTracerViewer::~G4RayTracerViewer() {}

void G4RayTracerViewer::SetView() {
  G4RayTracer* theTracer = 
    (G4RayTracer*) fSceneHandler.GetGraphicsSystem();
  if (fVP.GetFieldHalfAngle() == 0.) { // Orthogonal (parallel) projection.
    G4cout << 
      "G4RayTracerViewer::SetView: orthogonal projection not yet implemented."
	   << G4endl;
    return;
  }
  G4double bRadius =  // See G4ViewParameters for following procedure.
    fSceneHandler.GetScene() -> GetExtent().GetExtentRadius();
  if(bRadius<=0.) bRadius = 1.;
  G4double cDistance = fVP.GetCameraDistance(bRadius);
  G4ThreeVector eyePosition =
    fVP.GetCurrentTargetPoint() + cDistance * fVP.GetViewpointDirection();
  theTracer ->
    SetViewSpan(200. * fVP.GetFieldHalfAngle() / theTracer -> GetNColumn());
  theTracer -> SetTargetPosition(fVP.GetCurrentTargetPoint());
  theTracer -> SetEyePosition(eyePosition);
  theTracer -> SetHeadAngle(0.);  // Still to be figured out?????????
  theTracer -> SetLightDirection(-fVP.GetActualLightpointDirection());
}


void G4RayTracerViewer::ClearView() {}

void G4RayTracerViewer::DrawView() {
  SetView();  // Special graphics system - bypass ProcessView().
  G4RayTracer* theTracer = 
    (G4RayTracer*) fSceneHandler.GetGraphicsSystem();
  char fileName [100];
  G4std::ostrstream ost(fileName, 100);
  ost << "g4RayTracer." << fShortName << '_' << fFileCount++ << ".jpeg"
      << G4std::ends;
  theTracer -> Trace(fileName);
}
