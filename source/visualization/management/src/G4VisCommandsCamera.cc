// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsCamera.cc,v 1.5 2000-01-11 17:23:25 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "G4VisCommandsCamera.hh"

#include "G4UnitsTable.hh"

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
    G4cout << "Target point reset to centre of scene, ("
	   << G4BestUnit (pScene -> GetStandardTargetPoint ().x(), "Length")
	   << ", "
	   << G4BestUnit (pScene -> GetStandardTargetPoint ().y(), "Length")
	   << ", "
	   << G4BestUnit (pScene -> GetStandardTargetPoint ().z(), "Length")
	   << ")";
    G4cout << "\nZoom factor reset to 1.";
    G4cout << "\nDolly distance reset to 0.";
    G4cout << "\nViewpoint direction reset to +z.";
    G4cout << "\nUp vector reset to +y.";
    G4cout << G4endl;
    if (pVMan -> GetVerboseLevel () > 1) {
      pVMan -> PrintCurrentView ();
    }
  }
}
