// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VTreeViewer.cc,v 1.2 2001-06-15 07:12:37 stanaka Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4VTreeViewer.hh"

#include "G4ios.hh"
#include "g4std/strstream"

#include "G4VSceneHandler.hh"

G4VTreeViewer::G4VTreeViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name) {}

G4VTreeViewer::~G4VTreeViewer() {}

void G4VTreeViewer::SetView() {}

void G4VTreeViewer::ClearView() {}

void G4VTreeViewer::DrawView() {
  NeedKernelVisit ();  // Always need to visit G4 kernel.
  ProcessView ();
}
