// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML2FileViewer.hh,v 1.4 1999-12-15 14:54:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2FileViewer.hh
// Satoshi Tanaka and Yasuhide Sawada

#ifdef  G4VIS_BUILD_VRMLFILE_DRIVER

#ifndef G4VRML2FILE_VIEWER_HH
#define G4VRML2FILE_VIEWER_HH

#include "g4std/fstream"
#include "G4VViewer.hh"
#include "globals.hh"

class G4VRML2FileSceneHandler;

class G4VRML2FileViewer: public G4VViewer {
public:
	G4VRML2FileViewer(G4VRML2FileSceneHandler& scene, const G4String& name = "");
	virtual ~G4VRML2FileViewer();
	void ClearView();
	void DrawView();
	void ShowView();
	void FinishView();
private:
	void SetView(); // Do nothing. SendViewParameters will do its job.
	void SendViewParameters ()  ;

private:
	G4VRML2FileSceneHandler& fSceneHandler; // Reference to Graphics Scene for this view.
	G4std::ofstream&         fDest ;

	G4double      fViewHalfAngle ;	
	G4double      fsin_VHA       ;	

};

#endif //G4VRML2File_VIEW_HH
#endif //G4VIS_BUILD_VRMLFILE_DRIVER
