// Satoshi Tanaka  31th May 2001
// A dummy viewer for GAGTreeSceneHandler.

#include "G4GAGTreeViewer.hh"

#include "G4ios.hh"
#include "g4std/strstream"

G4GAGTreeViewer::G4GAGTreeViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  G4VTreeViewer(sceneHandler, name) {}

G4GAGTreeViewer::~G4GAGTreeViewer() {}
