// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1FileViewer.hh,v 1.4 1999-12-15 14:54:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML1FileViewer.hh
// Satoshi Tanaka and Yasuhide Sawada

#ifdef  G4VIS_BUILD_VRMLFILE_DRIVER

#ifndef G4VRML1FILE_VIEWER_HH
#define G4VRML1FILE_VIEWER_HH

#include "G4VViewer.hh"
#include "globals.hh"

class G4VRML1FileSceneHandler;

class G4VRML1FileViewer: public G4VViewer {
public:
	G4VRML1FileViewer(G4VRML1FileSceneHandler& scene, const G4String& name = "");
	virtual ~G4VRML1FileViewer();
	void ClearView();
	void DrawView();
	void ShowView();
	void FinishView();
private:
	void SetView();

private:
	G4VRML1FileSceneHandler& fSceneHandler; // Reference to Graphics Scene for this view.
};

#endif //G4VRML1File_VIEW_HH
#endif //G4VIS_BUILD_VRMLFILE_DRIVER
