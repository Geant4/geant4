// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML2FileView.hh,v 1.1 1999-01-07 16:15:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2FileView.hh
// Satoshi Tanaka and Yasuhide Sawada

#ifdef  G4VIS_BUILD_VRMLFILE_DRIVER

#ifndef G4VRML2FILE_VIEW_HH
#define G4VRML2FILE_VIEW_HH

#include <fstream.h>
#include "G4VView.hh"
#include "globals.hh"

class G4VRML2FileScene;

class G4VRML2FileView: public G4VView {
public:
	G4VRML2FileView(G4VRML2FileScene& scene, const G4String& name = "");
	~G4VRML2FileView();
	void ClearView();
	void DrawView();
	void ShowView();
	void FinishView();
private:
	void SetView(); // Do nothing. SendViewParameters will do its job.
	void SendViewParameters ()  ;

private:
	G4VRML2FileScene& fScene; // Reference to Graphics Scene for this view.
	ofstream&         fDest ;

	G4double      fViewHalfAngle ;	
	G4double      fsin_VHA       ;	

};

#endif //G4VRML2File_VIEW_HH
#endif //G4VIS_BUILD_VRMLFILE_DRIVER
