// Satoshi Tanaka  31th May 2001
// A dummy viewer for GAGTreeSceneHandler.

#ifndef G4GAGTREEVIEWER_HH
#define G4GAGTREEVIEWER_HH

#include "G4VTreeViewer.hh"

class G4GAGTreeViewer: public G4VTreeViewer {
public:
  G4GAGTreeViewer(G4VSceneHandler&,const G4String& name);
  virtual ~G4GAGTreeViewer();
};

#endif
