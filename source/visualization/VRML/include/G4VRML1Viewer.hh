// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1Viewer.hh,v 1.2 1999-01-11 00:48:00 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML1Viewer.hh
// Satoshi Tanaka & Yasuhide Sawada

#ifdef  G4VIS_BUILD_VRML_DRIVER

#ifndef G4VRML1_VIEWER_HH
#define G4VRML1_VIEWER_HH

#include "G4VViewer.hh"
#include "globals.hh"

class G4VRML1SceneHandler;

class G4VRML1Viewer: public G4VViewer {
public:
	G4VRML1Viewer(G4VRML1SceneHandler& scene, const G4String& name = "");
	~G4VRML1Viewer();
	void ClearView();
	void DrawView();
	void ShowView();
	void FinishView();
private:
	void SetView();

private:
	G4VRML1SceneHandler& fSceneHandler; // Reference to Graphics Scene for this view.
};

#endif //G4VRML1_VIEW_HH
#endif //G4VIS_BUILD_VRML_DRIVER
