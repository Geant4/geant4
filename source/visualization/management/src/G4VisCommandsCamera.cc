// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsCamera.cc,v 1.3 1999-01-11 00:48:31 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "G4VisCommandsCamera.hh"

void G4VisCommandCameraReset::SetValue () {
  G4VisManager* pVMan = G4VisManager::GetInstance ();
  if (pVMan -> IsValidView ()) {
    const G4Scene* pScene = pVMan -> GetCurrentScene ();
    G4ViewParameters& vp = pVMan -> SetCurrentViewParameters ();
    vp.SetCurrentTargetPoint (pScene -> GetStandardTargetPoint ());
    vp.SetZoomFactor (1.);
    vp.SetDolly (0.);
    vp.SetViewpointDirection (G4Vector3D (0., 0., 1.));
    vp.SetUpVector (G4Vector3D (0., 1., 0.));
    if (pVMan -> GetVerboseLevel () > 0) {
      G4cout << "Target point reset";
      G4cout << "\nZoom factor reset to 1.";
      G4cout << "\nDolly distance reset to 0.";
      G4cout << "\nViewpoint direction reset to +z.";
      G4cout << "\nUp vector set to +y.";
      G4cout << endl;
      if (pVMan -> GetVerboseLevel () > 1) {
	pVMan -> PrintCurrentView ();
      }
    }
  }
}
