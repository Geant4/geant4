// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1FileView.hh,v 1.1 1999-01-07 16:15:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML1FileView.hh
// Satoshi Tanaka and Yasuhide Sawada

#ifdef  G4VIS_BUILD_VRMLFILE_DRIVER

#ifndef G4VRML1FILE_VIEW_HH
#define G4VRML1FILE_VIEW_HH

#include "G4VView.hh"
#include "globals.hh"

class G4VRML1FileScene;

class G4VRML1FileView: public G4VView {
public:
	G4VRML1FileView(G4VRML1FileScene& scene, const G4String& name = "");
	~G4VRML1FileView();
	void ClearView();
	void DrawView();
	void ShowView();
	void FinishView();
private:
	void SetView();

private:
	G4VRML1FileScene& fScene; // Reference to Graphics Scene for this view.
};

#endif //G4VRML1File_VIEW_HH
#endif //G4VIS_BUILD_VRMLFILE_DRIVER
