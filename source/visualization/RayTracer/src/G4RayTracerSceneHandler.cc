// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayTracerSceneHandler.cc,v 1.1 2000-02-23 16:03:53 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4RayTracerSceneHandler.hh"

G4RayTracerSceneHandler::G4RayTracerSceneHandler(G4VGraphicsSystem& system,
						 const G4String& name):
  G4VSceneHandler(system, fSceneIdCount++, name)
{
  fSceneCount++;
}

G4RayTracerSceneHandler::~G4RayTracerSceneHandler()
{
  fSceneCount--;
}

G4int G4RayTracerSceneHandler::GetSceneCount() {
  return fSceneCount;
}

G4int G4RayTracerSceneHandler::fSceneIdCount = 0;

G4int G4RayTracerSceneHandler::fSceneCount = 0;
