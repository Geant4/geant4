// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayTracerViewer.hh,v 1.1 2000-02-23 16:03:52 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// John Allison  17th March 2000

#ifndef RAYTRACERVIEWER_HH
#define RAYTRACERVIEWER_HH

#include "G4VViewer.hh"

class G4RayTracerViewer: virtual public G4VViewer {
public:
  G4RayTracerViewer(G4VSceneHandler&,const G4String& name);
  virtual ~G4RayTracerViewer();
  void SetView();
  void ClearView();
  void DrawView();
private:
  G4int fFileCount;
};

#endif
