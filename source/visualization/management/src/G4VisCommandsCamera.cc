// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsCamera.cc,v 2.2 1998/07/23 02:18:55 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 

#include "G4VisCommandsCamera.hh"

void G4VisCommandCameraReset::SetValue () {
  G4VisManager* pVMan = G4VisManager::GetInstance ();
  if (pVMan -> IsValidView ()) {
    const G4SceneData& sd = pVMan -> GetCurrentSceneData ();
    G4ViewParameters& vp = pVMan -> SetCurrentViewParameters ();
    vp.SetCurrentTargetPoint (sd.GetStandardTargetPoint ());
    vp.SetZoomFactor (1.);
    vp.SetDolly (0.);
    if (pVMan -> GetVerboseLevel () > 0) {
      G4cout << "Target point reset";
      G4cout << "\nZoom factor reset to 1.";
      G4cout << "\nDolly distance reset to 0.";
      G4cout << endl;
      if (pVMan -> GetVerboseLevel () > 1) {
	pVMan -> PrintCurrentView ();
      }
    }
  }
}
