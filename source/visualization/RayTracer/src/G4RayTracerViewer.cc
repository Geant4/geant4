// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayTracerViewer.cc,v 1.4 2000-04-12 13:09:23 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4RayTracerViewer.hh"

#include "G4ios.hh"
#include "g4std/strstream"

#include "G4VSceneHandler.hh"
#include "G4Scene.hh"
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

  // Get radius of scene, etc.  (See G4OpenGLViewer::SetView().)
  // Note that this procedure properly takes into account zoom, dolly and pan.
  const G4Point3D& targetPoint = fVP.GetCurrentTargetPoint();
  G4double radius =  // See G4ViewParameters for following procedure.
    fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
  if(radius<=0.) radius = 1.;
  const G4double cameraDistance = fVP.GetCameraDistance(radius);
  const G4Point3D cameraPosition =
    targetPoint + cameraDistance * fVP.GetViewpointDirection().unit();
  const G4double nearDistance  = fVP.GetNearDistance(cameraDistance,radius);
  const G4double frontHalfHeight = fVP.GetFrontHalfHeight(nearDistance,radius);
  const G4double frontHalfAngle = atan(frontHalfHeight / nearDistance);

  // Calculate and set ray tracer parameters.
  theTracer->
    SetViewSpan(200. * frontHalfAngle / theTracer->GetNColumn());
  theTracer->SetTargetPosition(targetPoint);
  theTracer->SetEyePosition(cameraPosition);
  theTracer->SetHeadAngle(fVP.GetViewpointDirection().phi());
  theTracer->SetLightDirection(-fVP.GetActualLightpointDirection());
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
  theTracer->Trace(fileName);
}
