// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayXView.hh,v 2.1 1998/11/06 13:41:53 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Nikos Savvas  1st June 1997
// GEANT4 Ray Tracing X window view.

#ifdef G4VIS_BUILD_RAYX_DRIVER

#ifndef G4RAYXVIEW_HH
#define G4RAYXVIEW_HH

#include "G4RayView.hh"
#include "G4OwnColormapXView.hh"

class G4RayScene;

class G4RayXView: public G4RayView, public G4OwnColormapXView {
  
public:
  G4RayXView (G4RayScene& scene, const G4String& name = "");
  void DrawView ();

private:
  int ncells;
  /*****************************
  // An RGB map for now
  int red_max;
  int green_max;
  int blue_max;
  int red_mult;
  int green_mult;
  int blue_mult;
  ********************************/
};

#endif

#endif
